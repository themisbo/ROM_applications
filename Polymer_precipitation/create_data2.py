import numpy as np
import glob
from natsort import natsorted

# import matplotlib.pyplot as plt

all_paths = sorted(glob.glob("lmbda*"))

for path in range(len(all_paths)):
    all_tt = natsorted(glob.glob(str(all_paths[path]) + "/numpy_t*.npy"))
    arr_all = np.zeros((len(all_tt), 97, 97))
    for tt in range(len(all_tt)):
        test = np.load(str(all_paths[path]) + "/numpy_t" + str(tt + 1) + ".npy")
        test_dims = np.load(str(all_paths[path]) + "/numpy_dims" + str(tt + 1) + ".npy")
        test_order_inter = np.load(
            str(all_paths[path]) + "/numpy_order" + str(tt + 1) + ".npy"
        )
        test_order = test_order_inter[test_order_inter % 2 == 0]

        zz = test[test_order]
        xx = test_dims[:, 0][test_order_inter % 2 == 0]
        yy = test_dims[:, 1][test_order_inter % 2 == 0]

        # plt.scatter(xx, yy, c=zz, cmap=plt.get_cmap("coolwarm"))
        # plt.axis("scaled")

        arr_int1 = np.core.records.fromarrays([xx, yy, zz], names="x_val,y_val,conc")
        arr_int2 = np.unique(arr_int1, axis=0)
        arr_int3 = np.sort(arr_int2, order=["y_val", "x_val"])
        arr_fin = arr_int3["conc"]

        # plt.imshow(arr_fin.reshape(97, 97), cmap=plt.get_cmap("coolwarm"), origin="lower")
        arr_fin_sq = arr_fin.reshape(97, 97)
        arr_all[tt, :, :] = arr_fin_sq

    np.save(all_paths[path], arr_all)

# test = np.load("lmbda0.05.npy")
# for i in range(test.shape[0]):
#     plt.imshow(test[i], cmap=plt.get_cmap("coolwarm"), origin="lower")
#     plt.show()
