for i in 110 115 120 125 130 135 140 145; do
echo $i
combine -M Asymptotic ./muTau_SM0_mH$i.txt >> combine_0_$i.log
combine -M Asymptotic ./muTau_SM1_mH$i.txt >> combine_1_$i.log
combine -M Asymptotic ./muTau_SM2_mH$i.txt >> combine_2_$i.log
done
