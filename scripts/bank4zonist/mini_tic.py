

import numpy as np
from scipy.sparse import csr_matrix
from quantifs_methods import tic_normalization



# Example usage
if __name__ == "__main__":
    #
    intensity_matrix = np.array([[100.44, 400.55, 700.33],
                         [200.2,  500.7,  800.5 ],
                         [  0.  ,   0. ,    0.  ],
                         [300.4,  600.1,  900.3 ]])   # 3 features(peaks), 4 spots (each spectrum is a row (spot))

    intensity_matrix = csr_matrix(intensity_matrix)

    # Perform TIC normalization
    normalized_matrix = tic_normalization(intensity_matrix,
                                          ["s1", "s2", "s3", "s4"],
                                          [84.66, 88.31, 150.55])

    print("Original Intensity Matrix:")
    print(intensity_matrix.todense()) #, shape=(3,4)))
    print("\nTIC Normalized Intensity Matrix:")
    print(normalized_matrix)