set -euo pipefail


# Grid search for threeAFC models
for noise in 0 0.05 0.1 0.15 0.2 0.25 0.5 0.75 1.0 1.25 1.5
do 
    for encNoise in 0 0.25 0.5 0.75 1.0 1.25 1.5
    do
        for learningAdv in 0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25
        do 
            python multi_hub_rep.py --test threeAFC --jitter_pixels 2 --noise $noise --encoding_noise $encNoise --learning_adv $learningAdv
        done
    done
done

# Grid search for LE models
for noise in 0.25 0.5 0.75 1.0 1.25 1.5 
do 
    for learningAdv in 0 0.01 0.015 0.02 0.025 0.05 0.075 0.1
    do
        python multi_hub_rep.py --test learn_exemp --jitter_pixels 2 --noise $noise --learning_adv $learningAdv
    done
done

# Grid search for many_odd
for encNoise in 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.25 2.5 2.75 3.0 
do 
    python multi_hub_rep.py --test many_odd --encoding_noise $encNoise
done