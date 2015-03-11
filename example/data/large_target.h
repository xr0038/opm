namespace opm_test {

  std::vector<opm::xym> large_target = {
    { -7.76327e+02, -6.75951e+03, 1.47197e+03 },
    { 9.66459e+02, 4.65259e+03, 3.10357e+03 },
    { -9.09684e+02, 1.09579e+03, 8.89826e+02 },
    { -1.77178e+02, 2.72788e+03, 3.69029e+03 },
    { -8.59022e+02, -3.96860e+03, 2.60929e+03 },
    { 3.59723e+03, 1.99416e+02, 3.22774e+02 },
    { -3.18086e+02, 5.77775e+02, 1.34344e+03 },
    { 4.59024e+01, -4.10599e+03, 8.24477e+02 },
    { -8.08220e+01, -5.65336e+03, 2.07481e+03 },
    { 1.45325e+03, -4.44043e+03, 1.07851e+03 },
    { -1.81088e+03, 1.69310e+03, 1.89518e+03 },
    { 4.65854e+03, -7.59540e+02, 5.28705e+02 },
    { -3.46045e+03, 4.90306e+00, 2.27563e+03 },
    { 4.29489e+03, 1.76302e+03, 6.14934e+03 },
    { 5.17394e+03, -2.37197e+02, 5.93082e+02 },
    { -1.44043e+03, 5.19739e+03, 6.96992e+02 },
    { 4.29662e+03, 3.24240e+03, 1.40787e+03 },
    { 5.93404e+03, 8.45216e+02, 1.34379e+02 },
    { 3.18073e+03, 4.09337e+03, 1.50075e+03 },
    { -5.85575e+02, 4.83270e+03, 7.70423e+02 },
    { 1.80190e+03, 7.82615e+02, 1.73961e+03 },
    { 1.61187e+03, -1.33548e+02, 1.77832e+03 },
    { 6.30388e+02, -4.04372e+03, 5.37614e+02 },
    { -3.97349e+02, -1.18209e+03, 4.55592e+03 },
    { 4.14350e+03, -1.14969e+03, 2.60562e+02 },
    { 5.16646e+03, -9.05171e+02, 2.57712e+02 },
    { 2.23213e+03, 2.73605e+03, 1.59900e+03 },
    { 3.21927e+03, 2.37244e+03, 1.73091e+03 },
    { -1.29993e+03, -4.43575e+03, 6.92449e+02 },
    { -3.01718e+03, -4.07806e+03, 1.59411e+03 },
    { 4.51489e+03, 8.43359e+02, 2.16049e+03 },
    { -3.75913e+03, 2.69404e+03, 2.20387e+03 },
    { 4.25873e+03, -5.27197e+02, 8.36126e+02 },
    { -2.71431e+02, -7.85770e+02, 8.70143e+02 },
    { -3.64167e+03, -6.12434e+01, 1.14604e+03 },
    { -3.63822e+03, 2.21361e+03, 6.74100e+02 },
    { -3.73973e+03, -2.53561e+03, 6.43134e+02 },
    { -1.36385e+03, 3.51336e+03, 5.87109e+02 },
    { -2.56786e+03, 2.79511e+03, 2.29893e+03 },
    { -2.18200e+03, -5.12740e+03, 2.09539e+03 },
    { 2.43126e+03, 3.90420e+03, 2.24828e+03 },
    { 2.07263e+03, 2.53189e+03, 7.51483e+02 },
    { 1.64156e+03, 1.49270e+03, 3.31084e+03 },
    { 3.50110e+02, -9.93463e+02, 2.57179e+02 },
    { 6.65449e+02, -2.95836e+03, 4.50872e+01 },
    { -2.87562e+02, -1.17949e+03, 6.62756e+02 },
    { 7.61148e+02, -2.13928e+02, 2.09154e+03 },
    { -7.03237e+02, 4.17153e+03, 1.32938e+02 },
    { -2.85327e+03, 2.53431e+03, 2.62340e+03 },
    { 1.14160e+03, 4.47129e+03, 3.58916e+02 },
    { -1.81921e+03, -2.43824e+03, 2.32577e+03 },
    { -2.03463e+03, 3.49807e+03, 3.07569e+03 },
    { -4.64110e+03, 1.18325e+02, 1.28134e+02 },
    { 2.96371e+03, 1.37320e+03, 3.10353e+02 },
    { -3.22998e+03, 2.03976e+01, 1.08558e+03 },
    { -2.00991e+03, -2.47491e+03, 3.60371e+02 },
    { 1.78557e+03, 2.01374e+03, 2.43587e+03 },
    { -1.69231e+03, -2.78811e+03, 9.23898e+02 },
    { 7.33547e+02, -4.49193e+03, 9.09277e+02 },
    { -1.37733e+03, 4.46556e+03, 2.06546e+02 },
    { 6.27979e+03, 3.03776e+02, 1.58867e+03 },
    { 3.39511e+03, -1.07607e+03, 4.36506e+02 },
    { -3.34912e+03, 3.08809e+03, 1.37470e+03 },
    { 2.62067e+02, -3.53316e+03, 1.49378e+03 },
    { -3.22285e+03, 1.68533e+03, 7.95402e+02 },
    { 4.91017e+03, -9.96501e+02, 9.77110e+02 },
    { -2.20151e+03, -2.29305e+03, 1.51389e+03 },
    { 9.38684e+02, 2.89355e+03, 2.79614e+02 },
    { 1.87573e+03, -2.30467e+03, 9.31518e+02 },
    { -4.69143e+03, -2.22472e+03, 2.19962e+03 },
    { 1.20750e+03, -1.52006e+03, 1.79137e+03 },
    { -1.99386e+03, -5.37603e+02, 2.57534e+02 },
    { -4.23947e+03, -1.49506e+03, 9.49347e+02 },
    { 9.37009e+02, 4.73538e+03, 6.78562e+02 },
    { -5.13382e+03, 1.64842e+02, 4.77803e+03 },
    { 3.42251e+03, 2.22575e+03, 8.38735e+02 },
    { 5.62959e+02, 3.96661e+03, 1.28405e+03 },
    { -3.27590e+03, 5.02656e+01, 1.43661e+02 },
    { -8.23165e+02, 3.05357e+03, 8.40282e+02 },
    { 1.33377e+03, -3.31468e+03, 1.11781e+03 },
    { -1.54924e+03, 1.50533e+03, 7.62857e+02 },
    { 2.53177e+03, -3.14898e+03, 3.38085e+03 },
    { 1.58667e+03, -1.00961e+03, 8.50444e+02 },
    { -4.08742e+03, 2.81321e+03, 1.73948e+03 },
    { -3.80299e+03, 9.63898e+02, 1.17335e+03 },
    { 9.45706e+02, 2.72327e+03, 1.01096e+03 },
    { -6.22938e+02, -3.06383e+03, 4.22900e+02 },
    { -2.00880e+03, 1.83526e+03, 1.68424e+03 },
    { 4.26373e+03, 1.35848e+03, 2.57447e+03 },
    { 4.02017e+03, -6.47808e+02, 1.93066e+03 },
    { 1.02536e+03, 4.34384e+03, 2.01977e+03 },
    { 5.48370e+01, 3.93819e+03, 1.13833e+03 },
    { 2.44696e+03, -1.48650e+03, 7.89721e+02 },
    { 1.87718e+03, 5.94400e+03, 9.55162e+02 },
    { -2.05031e+02, 4.93508e+02, 1.72311e+02 },
    { -4.99055e+03, -1.78086e+03, 3.20092e+02 },
    { 3.29652e+03, 1.21053e+03, 1.58896e+03 },
    { 5.16467e+03, 1.49976e+03, 6.22425e+02 },
    { 1.51781e+03, 5.36700e+03, 1.31489e+03 },
    { -4.08836e+03, -2.19650e+03, 1.23724e+02 },
    { 4.58941e+03, 3.08375e+03, 5.44936e+02 },
    { -4.78038e+03, 1.61443e+03, 9.16959e+02 },
    { -4.90446e+03, -1.24041e+03, 2.28649e+03 },
    { 1.93128e+03, 1.09354e+03, 4.47444e+02 },
    { 3.50449e+02, 2.14649e+03, 5.43173e+02 },
    { 4.71809e+02, -1.82443e+03, 4.89667e+02 },
    { -4.85182e+03, 1.82466e+03, 1.72165e+03 },
    { -3.36143e+03, 1.82834e+03, 1.01058e+03 },
    { -3.37960e+03, 2.68858e+03, 1.69806e+03 },
    { -1.98166e+03, -3.85925e+03, 4.22452e+02 },
    { 1.78699e+03, -1.80045e+03, 1.16975e+03 },
    { -1.29849e+03, -6.01545e+03, 6.15806e+02 },
    { -1.71333e+03, 4.88871e+03, 2.72091e+03 },
    { -1.85713e+03, -4.05468e+03, 2.07736e+02 },
    { -1.59359e+03, -2.93700e+03, 6.95019e+02 },
    { 3.22970e+03, -2.95649e+03, 4.91881e+02 },
    { -2.18399e+03, -4.24615e+03, 8.45693e+02 },
    { 3.04760e+03, 1.72278e+03, 4.32420e+02 },
    { -2.30908e+03, 2.61587e+03, 9.82821e+02 },
    { 4.80179e+02, 5.83804e+03, 8.30081e+02 },
    { -6.09735e+03, 8.28435e+02, 2.63522e+02 },
    { 2.38656e+03, 2.60574e+03, 8.42756e+02 },
    { 2.24332e+03, -3.13835e+03, 7.87622e+02 },
    { 3.39027e+03, -1.70878e+03, 5.66544e+02 },
    { -8.81568e+02, -2.75432e+03, 6.41858e+02 },
    { -1.89347e+03, 4.29325e+03, 1.61396e+01 },
    { -4.24208e+03, 8.60299e+00, 5.44817e+02 },
    { 2.15036e+02, -5.73607e+03, 2.26100e+03 },
    { -2.28679e+03, 2.49775e+03, 1.37805e+03 },
    { -1.02167e+03, -5.17456e+03, 5.90629e+02 },
    { -4.19480e+02, -3.78029e+02, 9.21354e+02 },
    { -7.20443e+01, -3.39575e+03, 2.15543e+02 },
    { 1.64911e+03, 4.04495e+03, 5.20951e+02 },
    { -1.79584e+03, -2.81120e+03, 1.14702e+03 },
    { 5.29996e+03, 1.07634e+03, 1.08681e+02 },
    { -2.74471e+02, 5.75444e+03, 5.09307e+02 },
    { 2.96521e+03, -2.99089e+03, 6.14429e+02 },
    { -1.33254e+03, 1.44414e+02, 5.00832e+02 },
    { -5.36084e+02, -3.82964e+03, 8.24264e+02 },
    { 2.37721e+03, 2.87685e+03, 1.42534e+03 },
    { 4.10209e+03, 2.92306e+03, 5.09413e+02 },
    { -3.01505e+03, -1.45004e+03, 4.20762e+02 },
    { -4.67184e+03, -1.93333e+03, 5.94998e+01 },
    { 1.22545e+02, 1.07830e+03, 1.05084e+03 },
    { 3.25610e+02, -3.72206e+03, 1.54147e+03 },
    { -2.77524e+03, 4.51379e+02, 4.24534e+03 },
    { 4.60568e+03, 2.24685e+02, 1.08719e+02 },
    { 3.52433e+03, -7.50304e+02, 4.92978e+02 },
    { 3.42768e+03, 2.57315e+03, 1.18478e+03 },
    { 6.11467e+03, 7.76068e+02, 1.65584e+03 },
    { -5.51794e+02, 1.97222e+03, 2.80608e+03 },
    { 2.31985e+03, 5.69416e+03, 4.69608e+02 },
    { -3.74941e+03, 2.84976e+03, 2.20868e+02 },
    { -2.12298e+03, -5.25096e+03, 2.02019e+03 },
    { -1.43383e+03, 2.77638e+03, 1.24709e+02 },
    { 3.88396e+03, -7.00431e+02, 1.62494e+03 },
    { -3.28861e+03, -1.36556e+03, 9.65277e+02 },
    { 1.22032e+03, -4.63534e+03, 4.79475e+02 },
    { 5.31197e+03, 1.19705e+03, 7.96734e+01 },
    { -3.76391e+03, 2.84834e+03, 6.26418e+02 },
    { -1.00994e+03, -6.80758e+03, 2.60010e+03 },
    { 5.75898e+02, -4.30145e+03, 1.56260e+02 },
    { 1.34871e+02, -1.93463e+03, 4.90940e+02 },
    { 2.53894e+03, -8.61161e+01, 1.61816e+03 },
    { -1.70271e+03, -2.19391e+03, 1.26774e+03 },
    { -3.95434e+03, -8.43839e+01, 3.53409e+02 },
    { 1.10477e+03, -1.07326e+03, 5.00785e+02 },
    { 3.89158e+03, 3.15987e+03, 6.08815e+02 },
    { -1.30174e+03, 2.75754e+03, 2.88864e+03 },
    { 1.50463e+03, 6.40290e+03, 2.29726e+03 },
    { 4.13848e+03, 3.18584e+02, 2.42157e+02 },
    { 9.73320e+02, 3.55665e+03, 1.43449e+03 },
    { 3.94270e+03, 3.85757e+03, 2.14565e+03 },
    { 7.78331e+02, 5.99772e+03, 1.11057e+03 },
    { 7.58572e+02, -2.12078e+03, 3.25840e+02 },
    { -5.61337e+03, -1.02710e+03, 8.77036e+02 },
    { -1.10403e+03, 1.43974e+03, 8.67114e+02 },
    { 1.27753e+02, -5.91606e+03, 7.80293e+02 },
    { -3.29366e+03, -3.12430e+03, 4.14281e+03 },
    { 3.26559e+03, 1.95993e+03, 3.83618e+02 },
    { 3.22543e+03, 7.28115e+02, 1.04000e+03 },
    { -3.54937e+03, -1.20478e+03, 8.83276e+02 },
    { 9.18663e+02, 1.35395e+03, 1.05003e+03 },
    { 6.87707e+03, 6.20247e+02, 5.56984e+02 },
    { 1.52577e+03, 2.81254e+02, 2.52351e+03 },
    { -2.98094e+03, -8.95642e+02, 5.82003e+02 },
    { 2.97544e+03, -1.12286e+02, 3.39885e+03 },
    { -3.26717e+03, 4.24245e+02, 2.95361e+02 },
    { 1.80430e+03, 4.76727e+03, 9.11370e+02 },
    { -4.39466e+03, -1.94709e+03, 3.36484e+02 },
    { -5.69775e+03, 8.17545e+02, 7.90668e+01 },
    { 4.96628e+03, -1.59648e+02, 1.02907e+03 },
    { 2.16068e+03, -3.96570e+03, 6.44524e+02 },
    { 1.60950e+03, -2.29610e+03, 1.30478e+03 },
    { -5.29527e+03, 1.71062e+03, 2.22007e+03 },
    { -4.13016e+03, -6.72857e+02, 5.52975e+02 },
    { -2.70292e+03, -1.74287e+03, 3.24787e+02 },
    { -4.00121e+02, -5.04618e+03, 3.18392e+02 },
    { -3.98114e+03, 8.98176e+02, 1.52391e+03 },
    { -3.27351e+03, -1.80077e+03, 1.28805e+03 },
    { -3.80324e+01, 4.60334e+02, 5.36357e+02 },
    { 3.19376e+02, -5.32753e+03, 1.21905e+03 },
    { 3.08410e+03, 4.69499e+01, 2.89795e+03 },
    { -2.62367e+03, 1.81601e+02, 9.50709e+02 },
    { 1.35105e+03, 2.67368e+03, 6.04346e+02 },
    { -6.57123e+02, -6.60382e+03, 2.09033e+03 },
    { -6.46726e+03, -4.04816e+02, 7.08377e+02 },
    { -2.94270e+03, 2.42675e+03, 3.41742e+02 },
    { -1.77498e+03, -3.19508e+03, 9.49824e+02 },
    { -2.21031e+03, 4.82099e+02, 1.14757e+02 },
    { 2.35360e+02, 6.73255e+03, 1.22231e+03 },
    { 1.52278e+03, -3.63689e+03, 8.02423e+02 },
    { 1.75355e+03, -2.51517e+03, 1.13859e+03 },
    { -1.91521e+02, 3.13898e+03, 8.84145e+01 },
    { -1.19268e+03, -5.64915e+03, 5.14804e+02 },
    { -3.18378e+03, -3.96984e+02, 4.02227e+03 },
    { -5.39611e+01, -4.32135e+03, 7.46885e+02 },
    { -3.41867e+03, -1.98489e+03, 1.17453e+03 },
    { -2.84850e+03, -1.37608e+03, 2.20195e+03 },
    { -8.86442e+02, 8.18208e+02, 7.84710e+02 },
    { -6.14705e+01, 4.55514e+03, 3.09544e+03 },
    { -7.09468e+03, 1.07235e+01, 2.59609e+02 },
    { -7.99745e+02, -1.01789e+03, 1.18261e+03 },
    { 9.89780e+02, 3.65693e+03, 2.66772e+03 },
    { 2.60606e+03, 2.60918e+03, 1.17549e+03 },
    { -3.28741e+01, 5.49419e+03, 5.56784e+02 },
    { 1.23174e+03, 6.23130e+02, 2.13681e+02 },
    { -8.99279e+02, -5.32349e+03, 6.82168e+02 },
    { 3.27612e+03, 8.40928e+02, 7.17903e+02 },
    { -4.71629e+03, -2.64713e+03, 8.69026e+02 },
    { -4.05091e+03, -3.46092e+03, 1.16361e+03 },
    { -2.42534e+03, -1.51001e+03, 2.81188e+03 },
    { 1.19738e+03, 3.08145e+02, 9.17213e+02 },
    { -4.17737e+03, -3.38921e+03, 9.16227e+01 },
    { -2.96187e+03, -3.02223e+03, 1.68443e+03 },
    { -1.08989e+03, -5.81270e+03, 7.77459e+02 },
    { -1.75679e+03, 3.23345e+03, 1.98161e+03 },
    { 9.37642e+02, -3.65951e+03, 1.05032e+03 },
    { -3.91015e+03, -3.57614e+03, 4.51629e+02 },
    { 2.34106e+03, 4.71401e+03, 5.53794e+03 },
    { -2.28499e+03, 3.41771e+03, 4.90672e+03 },
    { -3.79475e+03, -1.03268e+03, 1.19368e+03 },
    { 3.37603e+03, 2.92176e+03, 4.97944e+02 },
    { -4.31376e+03, -2.62768e+03, 1.56186e+03 },
    { -1.61230e+03, -3.54788e+03, 4.23027e+03 },
    { 1.05096e+03, 6.32600e+03, 4.72549e+02 },
    { 2.69267e+03, -2.37991e+03, 1.45749e+03 },
    { 4.40505e+03, -6.43037e+02, 2.43711e+03 },
    { 1.07835e+03, 4.82492e+03, 9.68348e+02 },
    { -4.42712e+03, -6.83042e+02, 4.58086e+02 },
    { -3.95796e+02, 5.35901e+03, 7.67140e+02 },
    { -6.29329e+03, 7.82229e+02, 4.99693e+02 },
    { 2.48670e+03, -3.08421e+03, 8.24065e+02 },
    { 1.67657e+03, 7.62814e+02, 9.07687e+02 },
    { 2.72705e+03, -1.95082e+03, 2.67107e+03 },
    { 6.02200e+03, 1.52969e+03, 4.44963e+02 },
    { -4.40367e+02, -5.33023e+03, 3.77239e+02 },
    { -6.21844e+02, 1.54674e+03, 1.95277e+02 },
    { 3.74301e+02, 4.44863e+03, 2.20165e+03 },
    { -3.68777e+03, -2.22220e+03, 7.75177e+02 },
    { -4.91803e+02, 4.11137e+03, 2.40886e+03 },
    { -1.68749e+03, 4.64985e+02, 6.92710e+02 },
    { 3.85208e+03, 3.01234e+03, 2.45469e+03 },
    { 3.98161e+03, -1.24574e+03, 1.51343e+03 },
    { -8.70223e+02, -3.51507e+03, 7.80191e+02 },
    { 4.91702e+03, 2.88896e+03, 1.08752e+03 },
    { -2.91128e+03, -5.03777e+03, 8.34043e+02 },
    { -4.45297e+03, -2.25460e+03, 1.52165e+03 },
    { -2.18845e+03, -2.14825e+03, 9.27639e+02 },
    { -1.25428e+03, -4.84377e+03, 1.60338e+03 },
    { -9.75785e+01, -2.17328e+03, 9.19437e+02 },
    { -1.97381e+03, 3.20198e+03, 8.05955e+02 },
    { 4.08770e+03, -9.25001e+02, 6.79463e+02 },
    { -2.71130e+03, 2.91692e+03, 4.56956e+02 },
    { -1.50035e+03, -1.34470e+02, 2.03394e+03 },
    { 2.63594e+03, -3.32799e+03, 6.50154e+02 },
    { 5.29328e+02, -2.68460e+03, 1.94007e+03 },
    { -3.71435e+03, -5.17546e+01, 2.21197e+03 },
    { 1.73935e+03, 6.92218e+02, 1.17270e+03 },
    { -1.77430e+03, -9.30400e+02, 5.69441e+01 },
    { 2.38690e+03, 2.78255e+03, 2.95857e+02 },
    { 8.69851e+02, 2.64842e+03, 5.07729e+02 },
    { -2.13454e+03, -1.11156e+03, 2.18124e+03 },
    { 3.08565e+03, -1.21275e+03, 1.38655e+03 },
    { -2.57632e+03, 1.34692e+03, 1.21129e+03 },
    { -4.60591e+03, -7.76021e+02, 6.87900e+02 },
    { 2.06117e+03, 5.08972e+03, 1.83641e+03 },
    { -2.12815e+03, -5.40931e+03, 3.17246e+03 },
    { -2.74032e+03, -1.77280e+03, 2.56820e+03 },
    { 8.74163e+02, 2.25506e+03, 1.10056e+03 },
    { -2.88240e+03, -3.84480e+03, 1.04620e+03 },
    { 1.54890e+02, 2.48903e+03, 1.43082e+03 },
    { -1.64580e+03, 4.49854e+02, 2.01382e+03 },
    { -4.63602e+02, -4.30458e+03, 1.45350e+03 },
    { 2.53351e+03, 1.81485e+03, 2.04288e+03 },
    { -2.15613e+03, -2.04589e+03, 9.71574e+02 },
    { 3.43565e+03, 2.34590e+03, 9.10826e+02 },
    { -6.49830e+03, -6.12615e+01, 2.39954e+03 },
    { -7.61352e+02, -6.51970e+03, 3.18016e+02 },
    { -7.17043e+02, 5.08475e+03, 4.88542e+02 },
  };

}