import numpy as np

def conversion_matrix(old_t, new_t, p, tolerance):
    """
    """

    # Basis function matrix
    def R_matrix(x, mu, k, i, j):
        if (i > j) or (j > i+1):
            return 0
        
        if (i >=k) or (j > k):
            raise ValueError("indices out of range")
        if (j == i):
            return (old_t[mu+1+i] - x)/(old_t[mu+1+i]-old_t[mu+1+i-k])
        else :
            return (x - old_t[mu+1+i-k])/(old_t[mu+1+i]-old_t[mu+1+i-k])
        
    
    # init return value
    n_cps_old = old_t.size - p - 1
    n_cps_new = new_t.size - p - 1
    matrix = np.zeros((n_cps_old,n_cps_new))

    # Prepare
    following_n_equal_knots = np.zeros((new_t.size), dtype=np.int64)
    index_of_knot_span_base = np.zeros((new_t.size), dtype=np.int64)

    # Here we assume closed knot-vectors
    equal_index_count = 0
    current_old_id = old_t.size - 2
    current_knot_span = old_t.size - 1
    for i in range(new_t.size - 2, -1, -1):
        # Check if knots are equal
        if (np.abs(old_t[current_old_id] - new_t[i]) < tolerance):
            equal_index_count += 1
            if old_t[current_old_id] < old_t[current_old_id + 1]:
                current_knot_span = current_old_id
            current_old_id -= 1
        else : 
            if (old_t[current_old_id] > new_t[i]):
                raise ValueError("New knot spans is not subset")
            equal_index_count = 0

        # Assign
        if new_t[i] < old_t[current_knot_span]:
            index_of_knot_span_base[i] = current_knot_span - 1
        else:
            index_of_knot_span_base[i] = current_knot_span
        following_n_equal_knots[i] = equal_index_count
    # # print for debugging    
    print(following_n_equal_knots)
    print(index_of_knot_span_base)

    # 
    offset = 0
    for i in range(n_cps_old):
        mu = index_of_knot_span_base[i]

        # First non-zero matrix entry
        j = mu - p
        if (following_n_equal_knots[i+1] >= p):
            matrix[i, offset] = 1
            offset +=1
            continue
        # The p == 0 case is already caught with the previous exception

        matrix[i,j] = R_matrix(new_t[i + 1],mu,1,0,0)
        matrix[i,j+1] = R_matrix(new_t[i + 1],mu,1,0,1)
        
        # Loop over different matrices, last and first column of the d matrix
        # must be treated differently, because they only have one entry, we
        # loop backwards so we can avoid creating an additional auxiliary
        # vector.
        for d in range(2, p + 1):
            # First we assign the last column multiplication
            matrix[i,j+d] = R_matrix(new_t[i + d], mu, d, d, d-1) * matrix[i, j + d-1]
            # For matrix multiplication we can go backwards to avoid creating 
            # an additional aux vector
            for q in range(d-1, 0,-1):
                matrix[i,j + q] = R_matrix(new_t[i + d],mu,d,q-1,q) * matrix[i,j + q -1] + R_matrix(new_t[i + d],mu,d,q,q) * matrix[i,j + q]
    print(matrix)
# original_knot_vector = np.array([0,0,0,1,2,2,4,5,5,5])
# target_knot_vector = np.array([0,0,0,1,2,2,3,4,4,5,5,5])
# degree_original = 2
original_knot_vector = np.array([-1,-1,-1,-0,1,1,1])
target_knot_vector = np.array([-1,-1,-1,-1/2,0,1/2,1,1,1])
degree_original = 2

conversion_matrix(original_knot_vector, target_knot_vector, degree_original, 1e-12)
