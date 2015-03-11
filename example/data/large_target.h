namespace opm_test {

  std::vector<opm::xym> large_target = {
    { -3.77185e+02, -4.98694e+03, 1.46550e+03 },
    { 2.48674e+03, 1.17025e+03, 1.64598e+03 },
    { 7.06532e+02, 1.80401e+03, 9.88943e+01 },
    { 1.23335e+03, 5.00942e+03, 4.32086e+03 },
    { 5.61322e+03, -8.25497e+02, 1.86331e+03 },
    { 5.83091e+03, -1.21093e+02, 2.59116e+03 },
    { 8.39598e+02, 5.51971e+01, 2.87570e+03 },
    { 5.22478e+02, -2.74724e+03, 8.16361e+02 },
    { 2.88979e+03, 3.73722e+03, 4.33814e+02 },
    { 5.26074e+03, 1.51857e+03, 4.20057e+03 },
    { 5.33504e+03, -1.37489e+03, 8.70111e+02 },
    { 6.00412e+03, 4.25821e+03, 2.51777e+03 },
    { -1.73760e+03, -1.07056e+03, 2.07015e+03 },
    { -2.72320e+03, 1.81928e+03, 2.09370e+03 },
    { 3.24678e+03, 4.41132e+03, 3.75980e+03 },
    { 1.52194e+03, 1.60396e+03, 3.25187e+02 },
    { 1.56762e+03, 2.90207e+03, 1.36304e+02 },
    { 6.29826e+03, 2.69520e+03, 1.89759e+03 },
    { -4.79837e+03, 1.47387e+03, 7.67004e+02 },
    { -1.26134e+03, -3.76889e+03, 8.50636e+02 },
    { 2.74293e+03, 6.69404e+02, 7.71216e+02 },
    { 2.39315e+03, 3.12967e+03, 7.61851e+02 },
    { 1.92835e+03, 3.94229e+03, 6.21415e+02 },
    { 1.47533e+03, -1.12026e+03, 7.59936e+02 },
    { 1.63651e+03, 6.40290e+03, 3.91793e+03 },
    { -7.01709e+02, 1.47235e+03, 8.34839e+02 },
    { -2.44517e+03, 1.03647e+03, 6.28973e+02 },
    { 3.98598e+03, 1.11978e+03, 1.16860e+03 },
    { -2.72011e+03, 9.94440e+02, 2.17070e+03 },
    { 8.04946e+02, 4.93135e+03, 1.61882e+03 },
    { -1.17161e+03, 1.39782e+02, 5.98559e+02 },
    { -1.13229e+03, -1.03698e+03, 1.85155e+03 },
    { 4.10899e+03, 3.13104e+03, 4.65931e+02 },
    { -5.42642e+03, 4.27960e+02, 1.12732e+03 },
    { 6.18326e+03, -7.06637e+02, 1.80079e+03 },
    { 1.67382e+03, -5.06915e+03, 1.10568e+03 },
    { 2.45081e+03, 2.71147e+03, 2.66041e+01 },
    { 2.59480e+03, -1.77132e+03, 2.34779e+02 },
    { 2.01295e+03, 7.61961e+03, 1.31785e+03 },
    { 2.13781e+02, 5.20889e+03, 2.03791e+02 },
    { 6.91596e+03, 2.36287e+03, 6.57146e+02 },
    { 8.40493e+02, 4.65241e+03, 1.67195e+03 },
    { -4.03976e+02, -8.88371e+02, 8.95883e+02 },
    { 1.21233e+03, -5.53929e+03, 1.07320e+03 },
    { 5.80902e+03, -1.16100e+03, 5.35403e+02 },
    { 3.13864e+03, -2.28913e+03, 2.17162e+03 },
    { -2.25530e+03, -8.46890e+02, 9.54465e+01 },
    { 1.41457e+03, -1.21119e+03, 1.35176e+03 },
    { -8.57609e+02, -2.28879e+02, 1.65118e+02 },
    { 2.04633e+03, 3.86555e+03, 2.71310e+02 },
    { -3.12603e+02, -4.34872e+03, 6.95157e+02 },
    { 1.93223e+03, 3.49104e+03, 3.24356e+03 },
    { 6.15482e+03, 3.67358e+03, 4.54248e+02 },
    { 2.55849e+03, -1.30038e+03, 1.42613e+03 },
    { -2.33187e+03, 1.82808e+02, 1.06942e+03 },
    { 4.57657e+03, 1.41237e+03, 1.30421e+03 },
    { -9.60331e+02, 8.91834e+02, 7.09495e+02 },
    { 5.76664e+02, 6.31088e+03, 6.63298e+02 },
    { 3.68755e+02, 1.14552e+03, 4.97471e+01 },
    { -1.22903e+03, -9.81081e+02, 5.61602e+02 },
    { -4.21661e+02, 2.14054e+03, 1.19385e+03 },
    { 4.12142e+03, 3.75682e+03, 4.97618e+02 },
    { -8.71466e+02, -9.25183e+02, 1.34602e+03 },
    { 1.57163e+03, -2.60078e+03, 1.71735e+03 },
    { -6.76894e+02, 1.93224e+03, 9.32061e+02 },
    { 5.40538e+03, 3.00058e+03, 1.18727e+03 },
    { 2.55948e+03, -3.36769e+03, 2.31283e+03 },
    { 1.77790e+03, -2.45208e+03, 1.29786e+03 },
    { 4.94398e+03, -1.02638e+03, 1.32319e+03 },
    { 2.82078e+02, 3.88070e+03, 4.55203e+02 },
    { -2.14526e+02, 3.83921e+02, 1.12151e+03 },
    { 4.18256e+02, 1.87265e+03, 4.05198e+02 },
    { 2.05498e+03, 8.47000e+02, 8.39654e+02 },
    { -3.71107e+03, 1.12819e+03, 5.87335e+02 },
    { -4.27099e+03, -9.66991e+02, 1.57757e+03 },
    { -2.35262e+03, 1.48915e+03, 7.26293e+02 },
    { 1.37141e+03, -4.23489e+03, 6.25408e+02 },
    { 6.42729e+03, 3.41389e+03, 1.36998e+03 },
    { 6.07640e+02, -9.88490e+02, 1.19199e+03 },
    { -4.46067e+03, 1.30993e+03, 6.91148e+02 },
    { -1.53858e+03, 1.06299e+03, 2.07532e+03 },
    { -3.77981e+03, 2.92513e+02, 1.70628e+02 },
    { 5.15540e+03, 2.14090e+03, 1.06665e+03 },
    { 2.45755e+03, -3.87836e+03, 2.80362e+02 },
    { 8.22172e+03, 1.95912e+03, 1.03236e+03 },
    { 2.35798e+03, 4.20674e+03, 1.37808e+03 },
    { 4.26419e+01, 5.21394e+03, 1.36493e+03 },
    { 2.23817e+03, -3.39519e+03, 1.98154e+03 },
    { -8.74450e+02, 1.39565e+02, 5.89893e+02 },
    { -4.08706e+03, -1.13625e+03, 1.03157e+03 },
    { 3.94314e+03, -4.95490e+02, 8.79196e+02 },
    { 5.72190e+03, -6.08705e+02, 8.88384e+02 },
    { 5.21320e+03, 5.06377e+03, 1.48989e+03 },
    { 2.46843e+03, -1.84812e+02, 1.04244e+03 },
    { 1.10040e+03, 1.74844e+01, 8.89967e+01 },
    { 1.26894e+03, 3.47130e+03, 3.53978e+02 },
    { 2.70081e+03, 4.15011e+03, 2.59950e+02 },
    { 5.47742e+03, 2.36362e+03, 3.81211e+02 },
    { 3.72380e+03, 2.84507e+03, 2.61462e+02 },
    { -3.62032e+03, -5.60118e+02, 1.39789e+03 },
    { 1.47875e+03, 4.68965e+03, 1.64141e+03 },
    { 1.59002e+03, -4.34342e+03, 3.95805e+02 },
    { 3.78793e+02, 6.34826e+03, 9.21936e+02 },
    { -2.62579e+03, 2.87213e+03, 8.70739e+02 },
    { 9.33092e+02, 2.58490e+02, 5.93999e+02 },
    { 6.30006e+03, 1.66357e+03, 1.86771e+03 },
    { -2.64462e+03, 7.10588e+02, 1.37888e+03 },
    { -1.10143e+03, -3.06929e+03, 7.02208e+02 },
    { 3.28103e+03, -2.39294e+03, 1.43206e+03 },
    { 3.78277e+03, 4.25225e+03, 3.26163e+02 },
    { -2.48703e+03, -2.45885e+03, 5.71415e+02 },
    { -1.47121e+02, -5.47728e+03, 1.74073e+03 },
    { 3.37469e+03, 1.97076e+03, 1.24310e+03 },
    { -3.10805e+03, -2.54202e+03, 1.41138e+03 },
    { -3.37485e+03, 2.69477e+03, 2.24591e+02 },
    { -4.92740e+02, -4.21883e+03, 2.81793e+02 },
    { -2.93013e+03, 9.70186e+02, 1.55374e+03 },
    { -1.64290e+03, 9.96320e+02, 3.46410e+02 },
    { -5.02010e+03, 1.08098e+03, 8.30037e+02 },
    { -3.42401e+03, 2.55243e+03, 3.07355e+03 },
    { -8.52331e+02, -4.81232e+03, 5.47805e+02 },
    { 4.78382e+03, 1.66714e+03, 7.26687e+02 },
    { 6.59695e+03, 2.15887e+03, 2.68907e+03 },
    { -3.40284e+01, -2.59477e+03, 1.31460e+03 },
    { 1.64191e+02, -1.53133e+03, 1.66483e+02 },
    { 2.76077e+03, 2.76819e+02, 4.84392e+02 },
    { 1.52733e+03, -3.46249e+03, 1.42187e+03 },
    { 1.49288e+03, -4.58831e+03, 1.01013e+03 },
    { -2.87791e+03, 1.61842e+03, 2.88920e+02 },
    { 5.00931e+03, 4.32703e+03, 1.94944e+03 },
    { -3.49118e+03, 2.21919e+02, 1.04753e+03 },
    { -3.65628e+03, 2.01309e+03, 1.98720e+03 },
    { 5.12902e+03, 3.96548e+03, 9.61177e+02 },
    { 3.49885e+03, -5.89538e+02, 1.31931e+03 },
    { 9.24193e+02, 5.75435e+03, 2.16386e+03 },
    { 2.81103e+03, 7.19170e+03, 1.69087e+03 },
    { 1.79563e+02, 5.35015e+03, 1.65169e+03 },
    { -3.23378e+03, -7.29563e+02, 3.14218e+03 },
    { -2.68376e+03, -2.48084e+03, 4.72503e+02 },
    { -2.24861e+03, 2.82821e+03, 1.46013e+03 },
    { 1.09474e+03, 3.32165e+03, 1.14239e+03 },
    { -3.14170e+03, 8.96883e+02, 3.01965e+02 },
    { 4.66246e+03, 3.04239e+03, 4.42198e+02 },
    { 3.83179e+03, -1.57742e+03, 9.44736e+02 },
    { -2.53580e+02, -2.58090e+03, 4.45617e+02 },
    { 3.36106e+03, 3.53613e+03, 5.38066e+02 },
    { 2.01639e+03, -2.97533e+02, 7.60116e+02 },
    { -1.14740e+03, -2.72602e+03, 2.01727e+02 },
    { 1.17058e+02, -4.97558e+03, 2.95402e+02 },
    { 5.78439e+02, -2.73013e+03, 1.03288e+03 },
    { 2.00704e+03, 7.43676e+03, 1.46042e+02 },
    { 4.27708e+03, 5.19961e+03, 2.05160e+03 },
    { 1.40901e+03, 4.38026e+02, 1.35340e+03 },
    { 1.46878e+03, -2.91075e+03, 2.33004e+01 },
    { 5.51911e+03, -1.03634e+03, 2.60051e+03 },
    { -7.74561e+02, 1.36844e+03, 1.53568e+03 },
    { -8.80448e+02, -7.76876e+02, 2.70289e+03 },
    { 5.03069e+02, 2.10338e+03, 1.55493e+03 },
    { -1.89533e+02, 5.35246e+03, 7.32506e+02 },
    { 3.27815e+02, -4.72642e+03, 1.69222e+03 },
    { 4.17885e+03, -2.29804e+03, 1.95854e+02 },
    { 2.96994e+03, 5.59964e+03, 1.40217e+03 },
    { 2.70733e+03, 3.74897e+03, 1.55055e+02 },
    { 4.97633e+03, -2.01622e+03, 6.07486e+01 },
    { -2.43627e+03, -1.03377e+02, 5.13918e+02 },
    { -1.12013e+03, -1.11561e+03, 5.20876e+02 },
    { 3.27350e+03, -3.42332e+02, 1.30993e+03 },
    { 6.68673e+03, 1.67067e+03, 9.44218e+02 },
    { 5.59607e+03, 2.21295e+03, 9.73506e+02 },
    { 9.45674e+02, 5.26756e+03, 4.78678e+02 },
    { 2.99600e+01, 3.42119e+03, 1.62841e+03 },
    { 3.80756e+03, 5.91339e+03, 6.06222e+02 },
    { 3.09013e+03, -3.15520e+03, 1.59052e+03 },
    { -2.30603e+03, 2.37889e+03, 2.08328e+03 },
    { 1.96929e+03, 6.43186e+03, 4.78593e+02 },
    { 2.20373e+03, -1.53338e+03, 2.60177e+03 },
    { 1.92556e+03, -4.41317e+03, 1.84281e+03 },
    { 4.79680e+03, 1.13672e+03, 2.07304e+02 },
    { 3.07650e+03, -2.80590e+03, 1.67202e+03 },
    { -1.68023e+03, -3.66308e+03, 1.75179e+03 },
    { 5.37852e+03, 3.45736e+02, 4.46076e+03 },
    { 1.47833e+03, -3.38927e+03, 1.18116e+03 },
    { -3.42787e+03, -8.53977e+02, 2.01922e+02 },
    { -5.05889e+03, 2.72793e+02, 1.98488e+02 },
    { 4.48762e+03, 1.55601e+03, 6.63746e+02 },
    { 1.70497e+03, 8.30850e+02, 1.53968e+03 },
    { -2.67365e+02, 3.00940e+03, 2.60456e+02 },
    { 2.09296e+03, 4.72043e+03, 1.83924e+03 },
    { -1.44225e+03, -2.55784e+03, 5.87436e+03 },
    { 2.58613e+03, 7.72345e+03, 2.17023e+03 },
    { 3.46421e+03, 8.29261e+02, 1.74486e+03 },
    { 8.46959e+02, -1.99264e+03, 1.06387e+02 },
    { -2.36970e+03, -8.36444e+02, 1.62336e+03 },
    { -8.89279e+02, -3.56726e+03, 8.09844e+02 },
    { 3.16019e+02, -2.55076e+03, 1.16283e+03 },
    { -1.54763e+03, -3.31431e+03, 3.73704e+02 },
    { 2.28182e+03, -3.75961e+03, 7.23891e+02 },
    { 9.07016e+02, 6.33359e+03, 1.56419e+03 },
    { -1.99041e+03, -1.83250e+03, 1.90134e+03 },
    { 5.59762e+02, -2.00181e+03, 2.03761e+03 },
    { 2.29911e+03, -3.59366e+03, 4.62500e+02 },
    { -4.19363e+03, -1.36648e+03, 4.82039e+02 },
    { 2.53854e+03, -3.45625e+02, 3.81307e+02 },
    { 8.78407e+02, -3.18163e+03, 5.12961e+02 },
    { 3.48265e+03, 4.70408e+03, 5.80091e+02 },
    { 3.03295e+03, -8.05283e+02, 9.24923e+02 },
    { 1.53482e+00, 5.40395e+03, 2.71092e+03 },
    { 1.92959e+03, -7.71988e+02, 5.00320e+02 },
    { 7.46822e+02, 1.22057e+03, 2.72197e+02 },
    { 4.42797e+03, -2.32173e+03, 1.83804e+03 },
    { 1.92683e+03, -2.15790e+02, 1.04710e+03 },
    { 6.89325e+03, 3.01241e+03, 2.25484e+03 },
    { -1.87788e+03, -3.40474e+03, 3.78538e+02 },
    { 2.78985e+03, 6.09795e+03, 1.64358e+03 },
    { 5.08786e+03, 4.08099e+03, 6.27047e+02 },
    { 7.19272e+02, -2.77458e+03, 8.46149e+02 },
    { -2.16432e+03, -2.64326e+03, 7.02383e+02 },
    { 8.73387e+02, 2.58463e+03, 6.15728e+02 },
    { -3.79293e+02, 1.47955e+03, 4.10479e+03 },
    { -4.37003e+03, 1.64039e+03, 3.78236e+03 },
    { -2.55043e+03, 3.13837e+03, 4.05990e+02 },
    { 1.37891e+03, -3.71395e+03, 1.27653e+03 },
    { 1.81983e+03, -3.63498e+03, 2.76215e+03 },
    { 2.03144e+03, -3.39534e+03, 3.88470e+03 },
    { 1.23057e+03, 6.32073e+03, 4.78106e+02 },
    { 3.92701e+03, -2.34694e+03, 3.36222e+03 },
    { 1.56374e+03, 6.89351e+03, 2.31229e+03 },
    { 1.18814e+03, 2.04735e+03, 1.90509e+03 },
    { -1.07874e+03, -4.16630e+03, 3.28352e+03 },
    { 6.10842e+02, -1.10618e+03, 2.48630e+02 },
    { 5.47231e+03, 2.01676e+03, 8.34268e+02 },
    { 2.43065e+03, -3.80519e+03, 1.54693e+02 },
    { 3.71785e+03, -2.85664e+03, 3.23528e+03 },
    { 4.49668e+03, 1.60901e+03, 6.40558e+02 },
    { -4.21050e+02, 2.44644e+03, 3.49910e+03 },
    { 6.03161e+03, 1.70927e+03, 1.12734e+03 },
    { 3.25642e+03, 2.18908e+02, 9.12496e+01 },
    { 8.13463e+02, 3.89899e+01, 1.07799e+03 },
    { 4.02986e+01, -1.85439e+03, 9.03559e+02 },
    { 5.64934e+03, 9.19275e+02, 1.00288e+03 },
    { 3.98278e+03, -2.12181e+03, 1.40821e+03 },
    { 2.74481e+03, 4.70164e+03, 4.35639e+02 },
    { 4.30770e+03, -6.59531e+02, 1.42394e+03 },
    { 1.21823e+03, 5.49881e+03, 2.94097e+03 },
    { -1.85318e+03, -3.35223e+03, 1.92741e+02 },
    { 6.06445e+03, 2.89124e+03, 2.53350e+02 },
    { -3.67382e+02, -4.81517e+03, 3.38476e+02 },
    { -2.55197e+03, -2.14728e+03, 6.04552e+02 },
    { -3.33457e+01, 4.91382e+02, 1.33311e+03 },
    { 3.05414e+02, 6.15066e+03, 1.45258e+03 },
    { 3.42652e+03, -2.82561e+03, 8.43155e+02 },
    { -1.62920e+02, 4.41327e+03, 2.20464e+02 },
    { 5.36572e+03, -1.34998e+03, 9.14335e+02 },
    { 4.43308e+03, -8.26671e+02, 1.00992e+03 },
    { 6.68022e+03, 4.95519e+01, 3.18865e+02 },
    { 4.52923e+03, 1.33969e+03, 2.35591e+03 },
    { 2.09222e+03, 1.02034e+03, 4.40459e+02 },
    { 3.47993e+03, 2.63294e+03, 9.53018e+02 },
    { -3.15526e+03, -2.63932e+03, 2.52436e+03 },
    { -2.54775e+03, 3.61711e+03, 1.06506e+03 },
    { -2.48484e+03, -2.64087e+03, 1.74494e+03 },
    { 5.70707e+02, -4.72826e+03, 7.81436e+02 },
    { 3.30555e+01, 2.75071e+03, 2.43985e+03 },
    { -3.41539e+03, -7.03815e+02, 1.73997e+02 },
    { 5.56176e+02, -5.11650e+03, 4.85099e+02 },
    { -5.90624e+02, -5.17271e+03, 1.62865e+03 },
    { -2.98993e+03, -1.31928e+02, 7.84789e+02 },
    { 1.89629e+03, 5.07340e+03, 3.00230e+02 },
    { 1.04635e+03, 3.32289e+02, 6.31852e+02 },
    { 3.08503e+03, 3.75153e+03, 4.64957e+02 },
    { 4.13050e+03, -1.32664e+03, 3.27448e+02 },
    { 3.16176e+03, 3.89729e+02, 6.76347e+02 },
    { 2.55786e+03, 3.07867e+03, 1.33009e+03 },
    { 5.39139e+02, -1.12598e+03, 6.36561e+02 },
    { 4.23783e+03, 2.96646e+03, 1.80644e+03 },
    { 3.01828e+03, 3.34880e+03, 1.97737e+02 },
    { -1.18647e+03, -1.89384e+03, 1.04416e+03 },
    { 2.31880e+03, 4.48338e+03, 2.24167e+03 },
    { 2.56056e+03, -2.23829e+02, 2.16044e+03 },
    { -2.88117e+03, 3.19017e+03, 5.22225e+02 },
    { -7.61318e+02, -2.04059e+03, 4.58943e+02 },
    { -1.61634e+03, 1.83474e+03, 3.47493e+02 },
    { -4.70218e+01, -4.13165e+03, 3.86945e+02 },
    { 1.50812e+03, 2.59603e+03, 2.18530e+02 },
    { -6.50300e+02, -1.33245e+03, 4.50892e+02 },
    { 3.50571e+03, -4.64037e+02, 4.47823e+03 },
    { 2.90281e+02, 9.13549e+02, 1.26096e+02 },
    { 1.36253e+03, -4.93496e+03, 3.23540e+02 },
    { -7.43374e+02, -4.24213e+03, 7.83258e+02 },
    { 4.75408e+03, -8.60233e+02, 1.58411e+03 },
    { 6.52334e+03, 7.26121e+02, 1.85063e+03 },
    { 5.15884e+03, 2.50754e+03, 3.89886e+02 },
    { 5.78545e+03, 3.75123e+03, 1.15249e+03 },
    { -2.61725e+03, -2.21913e+03, 4.07490e+03 },
    { -9.97931e+02, -4.19217e+03, 1.36662e+03 },
    { 4.15993e+03, 2.28934e+03, 4.92461e+02 },
    { 6.98478e+03, 1.67949e+03, 9.90974e+02 },
    { 7.62009e+03, 1.34054e+03, 2.18514e+03 },
    { 4.54448e+03, -1.58630e+03, 5.40931e+02 },
    { 5.40659e+03, 2.74133e+03, 7.04866e+02 },
  };

}
