import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model, load_model

# import tensorflow.keras.backend as K
import sys, os
import glob
import pandas as pd
import numba as nb
import itertools
import json

from tensorflow.python.ops.gen_array_ops import parallel_concat


# sys.path.append("../python_scripts/svcca")
# import cca_core, pwcca

"""
Preprocessing functions
"""


@nb.jit(nopython=True)
def preprocess_rsaNumba(acts):
    assert acts.ndim in (2, 4)

    if acts.ndim == 4:
        imgs, x, y, chans = acts.shape
        newShape = x * y * chans
        newActs = acts.reshape((imgs, newShape))
        result = nb_cor(newActs, newActs)[0:imgs, imgs : imgs * 2]
    else:
        imgs = acts.shape[0]
        result = nb_cor(acts, acts)[0:imgs, imgs : imgs * 2]

    return result


@nb.jit(nopython=True)
def preprocess_peaRsaNumba(acts):
    # Alias of above to match others
    return preprocess_rsaNumba(acts)


@nb.jit(nopython=True, parallel=True)
def preprocess_speRsaNumba(acts):
    assert acts.ndim in (2, 4)

    if acts.ndim == 4:
        imgs, x, y, chans = acts.shape
        newShape = x * y * chans
        newActs = acts.reshape((imgs, newShape))
    else:
        imgs = acts.shape[0]
        newActs = acts

    result = nb_spearman(newActs, newActs)[0:imgs, imgs : imgs * 2]

    return result


# this
@nb.jit(nopython=True, parallel=True)
def preprocess_eucRsaNumba(acts):
    assert acts.ndim in (2, 4)

    if acts.ndim == 4:
        imgs, x, y, chans = acts.shape
        newShape = x * y * chans
        newActs = acts.reshape((imgs, newShape))
        # Preallocate array for RDM
        result = np.zeros((imgs, imgs))
        # Loop through each image and calculate euclidean distance
        for i in nb.prange(imgs):
            for j in nb.prange(imgs):
                result[i, j] = np.linalg.norm(newActs[i] - newActs[j])
                result[j, i] = result[i, j]
    else:
        imgs = acts.shape[0]
        # Preallocate array for RDM
        result = np.zeros((imgs, imgs))
        # Loop through each image and calculate euclidean distance
        for i in nb.prange(imgs):
            for j in nb.prange(imgs):
                result[i, j] = np.linalg.norm(acts[i] - acts[j])
                result[j, i] = result[i, j]

    return result


@nb.jit(nopython=True, parallel=True)
def nb_cov(x, y):
    # Concatenate x and y
    x = np.concatenate((x, y), axis=0)

    # Subtract feature mean from each feature
    for i in nb.prange(x.shape[0]):
        x[i, :] -= x[i, :].mean()

    # Dot product
    result = np.dot(x, x.T)

    # Normalization
    factor = x.shape[1] - 1
    result *= np.true_divide(1, factor)

    return result


@nb.jit(nopython=True, parallel=True)
def nb_cor(x, y):
    # Get covariance matrix
    c = nb_cov(x, y)

    # Get diagonal to normalize into correlation matrix
    d = np.sqrt(np.diag(c))

    # Divide by rows
    for i in nb.prange(d.shape[0]):
        c[i, :] /= d

    # Tranpose and divide it again
    c = c.T

    for i in nb.prange(d.shape[0]):
        c[i, :] /= d

    # Transpose back
    c = c.T

    return c


@nb.jit(nopython=True, parallel=True)
def nb_spearman(x, y):
    """
    Return Spearman rank correlation.
    """

    def _rank(x):
        idx = np.argsort(x)
        counter = 1
        for i in idx:
            x[i] = counter
            counter += 1

        return x

    # Get ranks for each vector
    xr = np.zeros(x.shape)
    yr = np.zeros(y.shape)
    for i in nb.prange(x.shape[0]):
        xr[i, :] = _rank(x[i, :])
        yr[i, :] = _rank(y[i, :])

    # Get rank correlation
    rankCor = nb_cor(xr, yr)

    return rankCor


# TODO: merge with interpolate
def preprocess_svcca(acts, interpolate=False):
    if len(acts.shape) > 2:
        acts = np.mean(acts, axis=(1, 2))
    # Transpose to get shape [neurons, datapoints]
    threshold = get_threshold(acts.T)
    # Mean subtract activations
    cacts = acts.T - np.mean(acts.T, axis=1, keepdims=True)
    # Perform SVD
    _, s, V = np.linalg.svd(cacts, full_matrices=False)

    svacts = np.dot(s[:threshold] * np.eye(threshold), V[:threshold])
    return svacts


# TODO: merge with interpolate
def preprocess_pwcca(acts, interpolate=False):
    if len(acts.shape) > 2:
        acts = np.mean(acts, axis=(1, 2))

    return acts


def get_threshold(acts):
    start = 0
    end = acts.shape[0]
    return_dict = {}
    ans = -1
    while start <= end:
        mid = (start + end) // 2
        # Move to right side if target is
        # greater.
        s = np.linalg.svd(
            acts - np.mean(acts, axis=1, keepdims=True), full_matrices=False
        )[1]
        # Note: normally comparing floating points is a bad bad but the precision we need is low enough
        if np.sum(s[:mid]) / np.sum(s) <= 0.99:
            start = mid + 1
        # Move left side.
        else:
            ans = mid
            end = mid - 1

    return ans


@nb.jit(nopython=True, parallel=True)
def preprocess_ckaNumba(acts):
    if acts.ndim == 4:
        nImg = acts.shape[0]
        nChannel = acts.shape[3]
        result = np.empty(shape=(nImg, nChannel), dtype="float32")
        for i in nb.prange(nImg):
            for j in nb.prange(nChannel):
                result[i, j] = np.mean(acts[i, :, :, j])

        return result
    else:
        return acts


"""
Correlation analysis functions
"""


# then this
@nb.jit(nopython=True, parallel=True)
def do_rsaNumba(rdm1, rdm2):
    """
    Pre: RDMs must be same shape
    """
    imgs = rdm1.shape[0]

    # Only use upper-triangular values
    upperTri = np.triu_indices(n=imgs, k=1)

    rdm1_flat = np.empty((1, upperTri[0].shape[0]), dtype="float32")
    rdm2_flat = np.empty((1, upperTri[0].shape[0]), dtype="float32")
    for n in nb.prange(upperTri[0].shape[0]):
        i = upperTri[0][n]
        j = upperTri[1][n]
        rdm1_flat[0, n] = rdm1[i, j]
        rdm2_flat[0, n] = rdm2[i, j]

    # Return pearson coefficient
    return nb_cor(rdm1_flat, rdm2_flat)[0, 1]


def do_svcca(acts1, acts2):
    """
    Pre: acts must be shape (neurons, datapoints) and preprocessed with SVD
    """
    svcca_results = cca_core.get_cca_similarity(
        acts1, acts2, epsilon=1e-10, verbose=False
    )
    return np.mean(svcca_results["cca_coef1"])


def do_pwcca(acts1, acts2):
    """
    Pre: acts must be shape (neurons, datapoints)
    """
    try:
        # acts1.shape cannot be bigger than acts2.shape for pwcca
        if acts1.shape <= acts2.shape:
            result = np.mean(pwcca.compute_pwcca(acts1.T, acts2.T, epsilon=1e-10)[0])
        else:
            result = np.mean(pwcca.compute_pwcca(acts2.T, acts1.T, epsilon=1e-10)[0])
    except np.linalg.LinAlgError as e:
        result = np.nan
        print(f"svd in pwcca failed, saving nan.")
        print(e)

    return result


@nb.jit(nopython=True)
def do_linearCKANumba(acts1, acts2):
    """
    Pre: acts must be shape (datapoints, neurons)
    """
    raise DeprecationWarning("Use do_linearCKANumba2 instead.")
    n = acts1.shape[0]
    centerMatrix = np.eye(n) - (np.ones((n, n)) / n)
    centerMatrix = centerMatrix.astype(nb.float32)

    # Top part
    centeredX = np.dot(np.dot(acts1, acts1.T), centerMatrix)
    centeredY = np.dot(np.dot(acts2, acts2.T), centerMatrix)
    top = np.trace(np.dot(centeredX, centeredY)) / ((n - 1) ** 2)

    # Bottom part
    botLeft = np.trace(np.dot(centeredX, centeredX)) / ((n - 1) ** 2)
    botRight = np.trace(np.dot(centeredY, centeredY)) / ((n - 1) ** 2)
    bot = (botLeft * botRight) ** (1 / 2)

    return top / bot


@nb.jit(nopython=True)
def do_linearCKANumba2(acts1, acts2):
    """
    Pre: acts must be shape (datapoints, neurons)
    """

    def _frobNorm(x):
        return np.sum(np.absolute(x) ** 2) ** (1 / 2)

    sim = _frobNorm(acts1.T @ acts2) ** 2
    normalization = _frobNorm(acts1.T @ acts1) * _frobNorm(acts2.T @ acts2)
    return sim / normalization


def correspondence_test(
    model1: str,
    model2: str,
    preproc_fun,
    sim_fun,
    names: str = None,
    rep_dir="../outputs/masterOutput/representations",
):
    """
    Return the results of the correspondence test given two model names using
    pregenerated representations in rep_dir with the preproc_fun and sim_funs.
    """
    # Get similarity function names
    if names is None:
        names = [fun.__name__ for fun in sim_fun]

    # Make directories for representations and get representations
    model1Glob = os.path.join(rep_dir, model1, f"{model1}l*.npy")
    model2Glob = os.path.join(rep_dir, model2, f"{model2}l*.npy")
    model1Glob = glob.glob(model1Glob)
    model2Glob = glob.glob(model2Glob)
    model1Glob.sort(key=lambda x: int(x.split("l")[-1].split(".")[0]))
    model2Glob.sort(key=lambda x: int(x.split("l")[-1].split(".")[0]))

    # Load representations
    model1Reps = [np.load(rep) for rep in model1Glob]
    model2Reps = [np.load(rep) for rep in model2Glob]

    # Get all representation layer combos and their similarities
    comboSims = list(itertools.product(range(len(model1Reps)), range(len(model1Reps))))
    comboSims = np.array(comboSims)
    comboSims = np.hstack((comboSims, np.zeros((len(comboSims), len(sim_fun)))))

    print("Generating representation simliarities", flush=True)
    for combo in comboSims:
        print(
            f"Finding the similarities between model 1 rep {int(combo[0])} and model 2 rep {int(combo[1])}",
            flush=True,
        )
        rep1 = model1Reps[int(combo[0])]
        rep2 = model2Reps[int(combo[1])]

        # Do analysis and save it
        output = multi_analysis(rep1, rep2, preproc_fun, sim_fun, names=names)
        combo[2:] = [output[fun] for fun in names]

    # Find the winner for each layer for each
    print("Finding winners for each layer", flush=True)
    winners = np.zeros((len(model1Reps), len(names)), dtype="int")
    for layer in range(len(model1Reps)):
        print(f"Finding the winner for layer {layer}")
        winners[layer, :] = np.argmax(comboSims[comboSims[:, 0] == layer, 2:], axis=0)

    winners

    return winners


def make_allout_model(model):
    """
    Creates a model with outputs at every layer that is not dropout.
    """
    inp = model.input

    modelOuts = [layer.output for layer in model.layers if "dropout" not in layer.name]

    return Model(inputs=inp, outputs=modelOuts)


def get_trajectories(directory, file_str="*", file_name=None):
    """
    Return a dataframe of the validation performance for each epoch for each
    model log in the directory, matches for a file_str if given. Saves to file
    if file_name is not None.
    """
    # Create dataframe
    colNames = ["weightSeed", "shuffleSeed", "epoch", "valAcc", "log"]
    out = pd.DataFrame(columns=colNames)

    logList = glob.glob(os.path.join(directory, file_str))

    for log in logList:
        logName = log.split("/")[-1]
        with open(log, "r") as f:
            lines = f.readlines()

        lines = [line.strip() for line in lines]

        # Go next if not finished.
        if not "final validation acc" in lines[-1]:
            continue

        # Find seeds
        snapshotLine = [line for line in lines if "Snapshot" in line][0]
        snapshotLine = snapshotLine.split()
        weightIndex = snapshotLine.index("weight") + 1
        shuffleIndex = snapshotLine.index("shuffle") + 1
        weightSeed = snapshotLine[weightIndex]
        shuffleSeed = snapshotLine[shuffleIndex]

        # Check for duplication
        seedCombos = out.groupby(["weightSeed", "shuffleSeed"]).groups.keys()
        if (str(weightSeed), str(shuffleSeed)) in seedCombos:
            print(f"Duplicate model parameter found from {logName}")

        # Get validation accuracy
        valAccs = [
            float(line.split("-")[-1].split(":")[-1].strip())
            for line in lines
            if "val_accuracy" in line
        ]
        nRows = len(valAccs)

        # Add to dataframe
        tmp = pd.DataFrame(
            data=[
                [weightSeed] * nRows,
                [shuffleSeed] * nRows,
                range(nRows),
                valAccs,
                [logName] * nRows,
            ],
            index=colNames,
        ).T
        out = out.append(tmp, ignore_index=True)

    if file_name is not None:
        out.to_csv(file_name, index=False)

    return out


def get_model_from_args(args, return_model=True, modelType="seed"):
    # Get model
    if hasattr(args, "model_name") and args.model_name == "mobilenet":
        model = tf.keras.applications.MobileNetV3Small(input_shape=(224, 224, 3))
        model.compile(metrics=["top_k_categorical_accuracy"])
        print(f"Model loaded: MobileNetV3Small", flush=True)
        model.summary()
        return model, "mobilenet", "."
    elif hasattr(args, "model_dir"):
        # List models in model_dir
        modelList = glob.glob(os.path.join(args.model_dir, "*.pb"))
        # Sort model list
        modelList.sort()
        # Get mode
        modelPath = modelList[args.model_index]
        modelName = modelPath.split("/")[-1].split(".")[0]
        model = load_model(modelPath)
    elif args.model_index is not None:
        import pandas as pd

        modelIdx = args.model_index

        # Load csv and get model parameters
        modelSeeds = pd.read_csv(args.model_seeds)
        weightSeed = modelSeeds.loc[modelSeeds["index"] == modelIdx, "weight"].item()
        shuffleSeed = modelSeeds.loc[modelSeeds["index"] == modelIdx, "shuffle"].item()

        # Load main model
        modelName = f"w{weightSeed}s{shuffleSeed}.pb"
        modelPath = os.path.join(args.models_dir, modelName)
        model = load_model(modelPath) if return_model else None
    elif args.shuffle_seed is not None and args.weight_seed is not None:
        weightSeed = args.weight_seed
        shuffleSeed = args.shuffle_seed

        # Load main model
        modelName = f"w{weightSeed}s{shuffleSeed}.pb"
        modelPath = os.path.join(args.models_dir, modelName)
        model = load_model(modelPath) if return_model else None
        print(f"Model loaded: {modelName}", flush=True)

    return model, modelName, modelPath


def get_funcs(method="all"):
    if method == "all":
        preprocFuns = [
            preprocess_peaRsaNumba,
            preprocess_eucRsaNumba,
            preprocess_speRsaNumba,
            preprocess_svcca,
            preprocess_ckaNumba,
        ]
        simFuns = [
            do_rsaNumba,
            do_rsaNumba,
            do_rsaNumba,
            do_svcca,
            do_linearCKANumba2,
        ]
        analysisNames = ["peaRsa", "eucRsa", "speRsa", "cca", "cka"]
    elif method == "rsa":
        preprocFuns = [
            preprocess_peaRsaNumba,
            preprocess_eucRsaNumba,
            preprocess_speRsaNumba,
        ]
        simFuns = [
            do_rsaNumba,
            do_rsaNumba,
            do_rsaNumba,
        ]
        analysisNames = ["peaRsa", "eucRsa", "speRsa"]
    elif method == "cs":
        preprocFuns = [
            preprocess_svcca,
            preprocess_ckaNumba,
        ]
        simFuns = [
            do_svcca,
            do_linearCKANumba2,
        ]
        analysisNames = ["svcca", "cka"]
    elif method == "good":
        preprocFuns = [preprocess_eucRsaNumba, preprocess_ckaNumba]
        simFuns = [do_rsaNumba, do_linearCKANumba2]
        analysisNames = ["eucRsa", "cka"]
    else:
        methods = method.split("-")
        preprocFuns = []
        simFuns = []
        analysisNames = []
        for string in methods:
            if string == "peaRsa":
                preprocFuns.append(preprocess_peaRsaNumba)
                simFuns.append(do_rsaNumba)
                analysisNames.append("peaRsa")
            elif string == "eucRsa":
                preprocFuns.append(preprocess_eucRsaNumba)
                simFuns.append(do_rsaNumba)
                analysisNames.append("eucRsa")
            elif string == "speRsa":
                preprocFuns.append(preprocess_speRsaNumba)
                simFuns.append(do_rsaNumba)
                analysisNames.append("speRsa")
            elif string == "cca":
                preprocFuns.append(preprocess_svcca)
                simFuns.append(do_svcca)
                analysisNames.append("cca")
            elif string == "cka":
                preprocFuns.append(preprocess_ckaNumba)
                simFuns.append(do_linearCKANumba2)
                analysisNames.append("cka")

    return preprocFuns, simFuns, analysisNames


def _split_comma_str(string):
    return [bit.strip() for bit in string.split(",")]


"""
Large scale analysis functions
"""


def multi_analysis(rep1, rep2, preproc_fun, sim_fun, names=None, verbose=False):
    """
    Perform similarity analysis between rep1 and rep2 once for each method as
    indicated by first applying a preproc_fun then the sim_fun. preproc_fun
    should be a list of functions to run on the representations before applying
    the similarity function at the same index in the list sim_fun. Elements of
    preproc_fun can be None wherein the representations are not preprocessed
    before being passed to the paired similarity function.
    """
    assert len(preproc_fun) == len(sim_fun)

    # If names is not none, just use the function names
    if names is None:
        names = [fun.__name__ for fun in sim_fun]

    # Loop through each pair
    simDict = {}
    counter = 0
    for preproc, sim in zip(preproc_fun, sim_fun):
        if verbose:
            print(f"___Preprocessing with {preproc.__name__}")

        rep1Copy = rep1.copy()
        rep2Copy = rep2.copy()

        # Preprocess each set of representations
        rep1Preproc = preproc(rep1Copy)
        rep2Preproc = preproc(rep2Copy)

        # Get similarity between reps
        try:
            # Print what we're doing
            if verbose:
                print(f"___Similarity with {names[counter]}", flush=True)

            simDict[names[counter]] = sim(rep1Preproc, rep2Preproc)
        except Exception as e:
            simDict[names[counter]] = np.nan
            print(f"{names[counter]} produced an error, saving nan.")
            print(e)

        counter += 1

        # Clean up memory
        del rep1Preproc, rep2Preproc, rep1Copy, rep2Copy

    return simDict


def get_reps_from_all(modelDir, dataset, outputDir=None):
    """
    Save representations for each model at every layer in modelDir by passing
    dataset through it.
    """
    # Get list of models
    models = os.listdir(modelDir)

    # Loop through models
    for model in models:
        print(f"Working on model: {model}")
        # Create allout model
        modelPath = os.path.join(modelDir, model)
        outModel = make_allout_model(load_model(modelPath))

        # Check if representation folder exists, make if not
        if outputDir is None:
            repDir = f"../outputs/masterOutput/representations/{model[0:-3]}"
        else:
            repDir = os.path.join(outputDir, model[0:-3])
        if not os.path.exists(repDir):
            os.makedirs(repDir)

        # Check if already done
        layerRepFiles = glob.glob(os.path.join(repDir, model[0:-3] + "l*"))
        if len(layerRepFiles) == len(outModel.outputs):
            print(
                "Layer representation files already exists, skipping.",
                flush=True,
            )
        else:
            # Get reps
            reps = outModel.predict(dataset, batch_size=128)

            # Save each rep with respective layer names
            for i, rep in enumerate(reps):
                np.save(f"{repDir}/{model[0:-3]}l{i}.npy", rep)


def get_unstruct_model_sims(
    repDir,
    layers,
    preprocFuns,
    simFuns,
    simNames,
    simMatType,
    outputDir="../outputs/masterOutput/similarities/",
    noise=None,
):
    """
    Return pairwise similarity from all model representations in repDir for the
    given layers using preprocFuns and simFuns.
    """
    # Get list of models
    models = glob.glob(os.path.join(repDir, "model*"))
    # Strip model directories
    models = [model.split("/")[-1] for model in models]

    # Loop through layers and combinations
    for layer in layers:
        # Create dataframe for model similarities
        sims = pd.DataFrame(columns=["model1", "model2"] + simNames)
        for combo in itertools.combinations(models, 2):
            print("Comparing models:", combo[0], combo[1])
            # Get reps
            rep1 = np.load(
                os.path.join(repDir, combo[0], combo[0] + "l" + layer + ".npy")
            )
            rep2 = np.load(
                os.path.join(repDir, combo[1], combo[1] + "l" + layer + ".npy")
            )

            if noise is not None:
                rep1 = rep1 + np.random.normal(
                    0, noise * np.std(rep1), size=rep1.shape
                ).astype("float32")
                rep2 = rep2 + np.random.normal(
                    0, noise * np.std(rep2), size=rep2.shape
                ).astype("float32")

            # Get similarities
            # can replace with this/then this
            simDict = multi_analysis(
                rep1, rep2, preprocFuns, simFuns, names=simNames, verbose=False
            )
            print(simDict)

            # Add to dataframe
            sims.loc[len(sims)] = list(combo) + list(simDict.values())

        # Save
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
        sims.to_csv(os.path.join(outputDir, f"{simMatType}_layer{layer}.csv"))

    return sims


def get_seed_model_sims(modelSeeds, repDir, layer, preprocFun, simFun, noise=None):
    """
    Return similarity matrix across all models in repDir and from a specific
    layer index using preprocFun and simFun. Representations should be in their
    own directory and use the following format within: w0s0l0.npy. The value
    after l is the layer index.
    """
    # Load csv and get model parameters
    import pandas as pd

    modelSeeds = pd.read_csv(modelSeeds)
    # Get all models and generate zeros matrix
    nModels = len(modelSeeds)
    modelSims = np.zeros(shape=(nModels, nModels))

    # Loop through, taking care to only do a triangle (assuming symmetry)
    for i, iRow in modelSeeds.iterrows():
        model1 = f"w{iRow.weight}s{iRow.shuffle}"
        print(f"==Working on index {i} model: {model1}==")
        # Load representations of i and preprocess
        rep1 = np.load(os.path.join(repDir, model1, f"{model1}l{layer}.npy"))

        if noise is not None:
            rep1 = rep1 + np.random.normal(
                0, noise * np.std(rep1), size=rep1.shape
            ).astype("float32")
        rep1 = preprocFun(rep1)
        for j, jRow in modelSeeds.iterrows():
            if i > j:  # Only do triangle
                continue

            model2 = f"w{jRow.weight}s{jRow.shuffle}"
            print(f"-Calculating similarity against index {j} model: {model2}")

            # Load representations of j and preprocess
            rep2 = np.load(os.path.join(repDir, model2, f"{model2}l{layer}.npy"))
            if noise is not None:
                rep2 = rep2 + np.random.normal(
                    0, noise * np.std(rep2), size=rep2.shape
                ).astype("float32")
            rep2 = preprocFun(rep2)

            # Do similarity calculation
            tmpSim = simFun(rep1, rep2)
            modelSims[i, j] = tmpSim
            modelSims[j, i] = tmpSim

    return modelSims


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Perform some type of analysis, intended to be used in HPC"
    )
    parser.add_argument(
        "--analysis",
        "-a",
        type=str,
        help="type of analysis to run",
        choices=["correspondence", "getReps", "seedSimMat", "itemSimMat"],
    )
    parser.add_argument(
        "--model_index",
        "-i",
        type=int,
        help="model index to select weight and shuffle seeds",
    )
    parser.add_argument(
        "--shuffle_seed", type=int, help="shuffle seed of the main model"
    )
    parser.add_argument("--weight_seed", type=int, help="weight seed of the main model")
    parser.add_argument(
        "--model_seeds",
        type=str,
        default=None,
        help="file location for csv file with model seeds",
    )
    parser.add_argument(
        "--dataset_file",
        type=str,
        default="../outputs/masterOutput/dataset.npy",
        help="npy file path for the image dataset to use for analysis",
    )
    parser.add_argument(
        "--models_dir",
        type=str,
        default="../outputs/masterOutput/models",
        help="directory for all of the models",
    )
    parser.add_argument(
        "--reps_dir",
        type=str,
        default="../outputs/masterOutput/representations",
        help="directory for representations",
    )
    parser.add_argument(
        "--simSet",
        type=str,
        default="all",
        help="which set of similarity functions to use",
    )
    parser.add_argument(
        "--layer_index",
        "-l",
        type=_split_comma_str,
        help="which layer to use, must be positive here, split by comma",
    )
    parser.add_argument(
        "--simMatType",
        type=str,
        help="what to name the output similarity matrix file",
    )
    parser.add_argument("--output_dir", "-o", type=str, default=None)
    parser.add_argument(
        "--noise",
        type=float,
        default=None,
        help="if set, add this proportion of noise to the representations based on their standard deviation",
    )
    args = parser.parse_args()

    np.random.seed(2022)
    # Now do analysis
    if args.analysis == "correspondence":
        print("Performing correspondence analysis.", flush=True)

        preprocFuns, simFuns, analysisNames = get_funcs(args.simSet)

        # Get model
        if args.model_seeds is not None:
            _, modelName, _ = get_model_from_args(args, return_model=False)
            modelName = modelName.split(".")[0]
        else:
            modelList = glob.glob(os.path.join(args.reps_dir, "*.npy"))
            # Get mode
            modelPath = modelList[args.model_index]
            modelName = modelPath.split("/")[-1].split(".")[0].split("-Reps")[0]

        # List model representations and make combinations
        reps = glob.glob(args.reps_dir + "/*")
        for rep in reps:
            if "w" in rep and "s" in rep:
                reps = [rep.split("/")[-1] for rep in reps if "w" in rep and "s" in rep]
            else:
                reps = [
                    rep.split("/")[-1].split(".")[0].split("-Reps")[0] for rep in reps
                ]
        # reps = [
        #     rep.split("/")[-1] for rep in reps if "w" in rep and "s" in rep else
        # ]
        repCombos = list(itertools.combinations(reps, 2))
        repCombos = [x for x in repCombos if x[0] == modelName]

        # Prepare dataframes
        if args.output_dir is not None:
            fileName = f"{args.output_dir}/{modelName}Correspondence.csv"
        else:
            fileName = (
                f"../outputs/masterOutput/correspondence/{modelName}Correspondence.csv"
            )
        if os.path.exists(fileName):
            # Load existing dataframe
            print("path exists")
            winners = pd.read_csv(fileName, index_col=0)
        else:
            numLayers = len(glob.glob(f"{args.reps_dir}/{reps[0]}/{reps[0]}l*.npy"))
            winners = pd.DataFrame(
                sum([[combo] * numLayers for combo in repCombos], []),
                columns=["model1", "model2"],
            )
            winners[analysisNames] = -1

        # Find the winners
        for model1, model2 in repCombos:
            print(f"Comparing {model1} and {model2}", flush=True)
            if np.all(
                winners.loc[
                    (winners["model1"] == model1) & (winners["model2"] == model2),
                    analysisNames,
                ]
                == -1
            ):
                winners.loc[
                    (winners["model1"] == model1) & (winners["model2"] == model2),
                    analysisNames,
                ] = correspondence_test(
                    model1, model2, preprocFuns, simFuns, names=analysisNames
                )

                print("Saving results", flush=True)
                winners.to_csv(fileName)
            else:
                print("This pair is complete, skipping", flush=True)

    elif args.analysis == "getReps":
        print("Getting representations each non-dropout layer", flush=True)

        # Load dataset
        print("Loading dataset", flush=True)
        dataset = np.load(args.dataset_file)
        print(f"dataset shape: {dataset.shape}", flush=True)

        # Run it!
        get_reps_from_all(args.models_dir, dataset, args.output_dir)
    elif args.analysis == "seedSimMat":
        print("Creating model similarity matrix.", flush=True)
        preprocFuns, simFuns, simNames = get_funcs(args.simSet)

        for layer in args.layer_index:
            for preprocFun, simFun, simName in zip(preprocFuns, simFuns, simNames):
                print(f"Working on layer {layer} with {simFun.__name__}")
                simMat = get_seed_model_sims(
                    args.model_seeds,
                    args.reps_dir,
                    layer,
                    preprocFun,
                    simFun,
                    args.noise,
                )

                if args.output_dir is None:
                    np.save(
                        f"../outputs/masterOutput/similarities/simMat_l{layer}_{simName}.npy",
                        simMat,
                    )
                else:
                    if not os.path.exists(args.output_dir):
                        os.makedirs(args.output_dir)
                    np.save(
                        os.path.join(args.output_dir, f"simMat_l{layer}_{simName}.npy"),
                        simMat,
                    )
    elif args.analysis == "itemSimMat":
        print(
            "Creating model similarity matrix for item weighting differences.",
            flush=True,
        )
        preprocFun, simFun, simNames = get_funcs(args.simSet)

        get_unstruct_model_sims(
            args.reps_dir,
            args.layer_index,
            preprocFun,
            simFun,
            simNames,
            args.simMatType,
            args.output_dir,
            args.noise,
        )
    else:
        import time
        # Simulate pairs of representations
        nRepeats = 10
        n = 200
        rep1K = 512
        rep2K = 256

        # Simulate random representations
        rep1 = np.random.normal(0, 1, size=(n, rep1K)).astype(np.float32)
        rep2 = np.random.normal(0, 1, size=(n, rep2K)).astype(np.float32)

        # Calculate RDM
        def numbaRSA(rep1, rep2):
            rdm1 = preprocess_eucRsaNumba(rep1)
            rdm2 = preprocess_eucRsaNumba(rep2)

            # Calculate RSA
            print(do_rsaNumba(rdm1, rdm2))

        t0 = time.time()

        for _ in range(nRepeats):
            numbaRSA(rep1, rep2)

        print(f"Time taken: {time.time() - t0}")
