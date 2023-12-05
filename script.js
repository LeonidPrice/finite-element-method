
var matrix = {
    errors : [
        "MATRIX ERROR: The matrices aren't compatible.\nThe number of columns of matrix A must be equal to the number of rows of matrix B.",
        "MATRIX ERROR: The matrices aren't compatible.\nThe number of rows and columns of matrices A and B must be equal.",
        "MATRIX ERROR: Negative Degree.\nThe degree of the matrix must be positive and greater than one.",
        "MATRIX ERROR: Matrices are not equal.",
        "MATRIX ERROR: The system is undefined and has an infinite number of solutions.",
        "MATRIX ERROR: The system is incompatible and has no solutions.",
        "MATRIX ERROR: The number of iterations is not enough to obtain a solution."
    ],

    /**
     * Sorting an array of numbers.
     * 
     * O(n^2)
     * @param {Number[]} array Unsorted array.
     * @returns {Number[]} Sorted array.
     */
    selection_sort: function (array) {
        var sorted_array = [];
        var index = [];

        while (array.length > 0) {
            var largest = array[0];
            var largest_index = 0;

            for (var i = 1; i < array.length; i++) {
                if (array[i] > largest) {
                    largest = array[i];
                    largest_index = i;
                }
            }

            index.push(largest_index);
            sorted_array.push(array.splice(largest_index, 1)[0]);
        }

        return [sorted_array, index[0]];
    },

    /**
     * Creates a unit matrix nxn of the form:
     * 
     *     1  0  0  
     *     0  1  0    at n = 3
     *     0  0  1
     * 
     * O(n^2)
     * @param {Number} n Dimension of the matrix.
     * @returns {Number[][]} E.
     */
    identity: function(n) {
        var E = [];

        for (var i = 0; i < n; i++) {
            E[i] = [];
            for (var j = 0; j < n; j++) {
                E[i][j] = i === j ? 1 : 0;
            }
        }
        return E;
    },

    /**
     * Calculates the rank of the matrix.
     * 
     * O(min(m,n) * m * n) where:
     * 
     * n - count of rows;
     * 
     * m - count of columns.
     * @param {Number[][]} A Source Matrix.
     * @param {Number} epsilon Error magnitude. Default: 1e-6.
     * @returns {Number} Rank of the matrix.
     */
    rank: function(A, epsilon = 1e-6) {
        if (!A || A.length === 0 || A[0].length == 0) {
            return false;
        }

        var line_used = [];
        var rank = A.length;

        A.length > A[0].length || A.length == A[0].length ? rank = A.length : rank = A[0].length;

        for (var i = 0; i < A.length; i++) {
            line_used[i] = false;
        }
    
        for (var i = 0; i < A[0].length; i++) {
            var j;
            for (var j = 0; j < A.length ; j++) {
                if (!line_used[j] && Math.abs(A[j][i]) > epsilon) {
                    break;
                }
            }
            if (j == A.length) {
                rank--;
            } else {
                line_used[j] = true;
                for (var p = i + 1; p < A[0].length; p++) {
                    A[j][p] /= A[j][i];
                }
                for (var k = 0; k < A.length; k++) {
                    if (k != j && Math.abs(A[k][i]) > epsilon) {
                        for (var p = i + 1; p <= A[0].length - 1; p++) {
                            A[k][p] -= A[j][p] * A[k][i];
                        }
                    }
                }
            }
        }
        return rank;
    },

    /**
     * Computes the first and infinite norms of a matrix.
     * 
     * O(n^2)
     */
    norma: {
        /**
         * O(n^2)
         * @param {Number[][]} A Source matrix.
         * @returns {Number} First norm of the matrix.
         */
        first: function(A) {
            if (!A || A.length === 0 || A[0].length == 0) {
                return false;
            }

            var first_norma = [];
            var s = 0;

            for (var i = 0; i < A.length; i++) {
                for (var j = 0; j < A[0].length; j++) {
                    s += Math.abs(A[j][i]);
                    
                }
                first_norma.push(s)
                s = 0;
            }
            return matrix.selection_sort(first_norma)[0][0];
        },

        /**
         * O(n^2)
         * @param {Number[][]} A Source matrix.
         * @returns {Number} Infinity norm of the matrix.
         */
        infinity: function(A) {
            if (!A || A.length === 0 || A[0].length == 0) {
                return false;
            }

            var infinity_norma = [];
            var s = 0;

            for (var i = 0; i < A.length; i++) {
                for (var j = 0; j < A[0].length; j++) {
                    s += Math.abs(A[i][j]);
                    
                }
                infinity_norma.push(s)
                s = 0;
            }
            return matrix.selection_sort(infinity_norma)[0][0];
        }
    },

    /**
     * Addition compatible matrices.
     * 
     * O(n^2)
     * @param {Number[][]} A Matrix A.
     * @param {Number[][]} B Matrix B. 
     * @returns {Number[][]} Matrix A+B.
     */
    add: function(A, B) {
        if (!A || !B || A.length === 0 || B.length === 0 || A.length !== B.length || A[0].length !== B[0].length) {
            console.error(this.errors[1]);
            return false;
        } else {
            var AB = [];
            for (var i = 0; i < A.length; i++) {
                AB[i] = [];
                for (var j = 0; j < A[0].length; j++) {
                    AB[i][j] = A[i][j] + B[i][j];
                }
            }
            return AB;
        }
    },

    /**
     * Subtraction compatible matrices.
     * 
     * O(n^2)
     * @param {Number[][]} A Matrix A.
     * @param {Number[][]} B Matrix B. 
     * @returns {Number[][]} Matrix A-B.
     */
    subtract: function(A, B) {
        if (!A || !B || A.length === 0 || B.length === 0 || A.length !== B.length || A[0].length !== B[0].length) {
            console.error(this.errors[1]);
            return false;
        } else {
            var AB = [];
            for (var i = 0; i < A.length; i++) {
                AB[i] = [];
                for (var j = 0; j < A[0].length; j++) {
                    AB[i][j] = A[i][j] - B[i][j];
                }
            }
            return AB;
        }
    },

    /**
     * Multiplies a matrix by a number.
     * 
     *     5  1  8             10  2  16
     *     3  2  0   will be    6  4   0  at n = 2
     *     6  7  4             12  14  8
     * 
     * O(n^2)
     * @param {Number[][]} A Matrix A.
     * @param {Number} n Any number. 
     * @returns {Number[][]} A*n.
     */
    scalar: function(A, n) {
        var A_n = [];

        for (var i = 0; i < A.length; i++) {
            A_n[i] = [];
            for (var j = 0; j < A[0].length; j++) {
                A_n[i][j] = A[i][j] * n;
            }
        }
        return A_n;
    },

    /**
     * Elevates the matrix A to positive degree n.
     * 
     * O(n * i * j * k) where:
     * 
     * n - value of the degree;
     * 
     * i - amount of rows in A;
     * 
     * j - amount of columns in B;
     * 
     * k - amount of rows in B.
     * @param {Number[][]} A Matrix A. 
     * @param {Number} n Degree of the matrix.
     * @returns {Number[][]} A^n
     */
    degree: function(A, n) {
        if (n <= 0) {
            console.error(this.errors[2]);
            return;
        }
        if (n == 1) {
            return A;
        }

        var A_previous = this.multiply(A, A);

        if (n == 2) {
            return A_previous;
        } else {
            var A_n;
            for (var i = 0; i <= n - 3; i++) {
                A_n = this.multiply(A, A_previous);
                A_previous = A_n;
            }
            return A_n;
        }
    },

    /**
    * Solves the system of equations AX = B with respect to X. SOR method.
    * 
    * O(k * n^2) where:
    * 
    * k - amount of iterations;
    * 
    * n - size of the matrix A.
    * @param {Number[][]} A Matrix A.
    * @param {Number[][]} B Matrix B.
    * @param {Number} epsilon Error magnitude. Default: 1e-6.
    * @param {Number} omega Relaxation parameter in the range 0 ... 1. Default: 0.5.
    * @param {Number} max_iterations Maximum number of iterations. Default: 1000.
    * @returns {Number[][]} Vector of roots X.
    */
    roots: function (A, B, epsilon = 1e-6, omega = 0.5, max_iterations = 1000) {
        if (!A || !B || A.length === 0 || B.length === 0) {
            return false;
        }
        
        var X0 = [];
        var X_tau = []
        var X = [];
        var error = [];

        var k = 0;
        while (k <= max_iterations) {
            for (var i = 0; i < B.length; i++) {
                X0[i] = X[i] || [0];
                X_tau[i] = X[i] || [0];
                X[i] = [0];
            }
        
            for (var i = 0; i < B.length; i++) {
                var S = 0;
                for (var j = 0; j < A.length; j++) {
                    if (i !== j) {
                        S += A[i][j] * X_tau[j]
                    }
                }
    
                X_tau[i] = 1 / A[i][i] * (B[i] - S);
    
                for (var j = 0; j < A.length; j++) {
                    X[j][0] = omega * X_tau[j] + (1 - omega) * X0[j][0];
                }
    
                error.push(Math.abs(X[i] - X0[i]))
            }
            
            error = this.selection_sort(error)[0][B.length - 1];
    
            if (error < epsilon) {
                return X;
            }
            error = [];
            k++;
        }
        console.error(this.errors[6]);
        return false;
    },

    /**
     * Transpose matrix A.
     * 
     * O(i * j) where:
     * 
     * i - amount of rows;
     * 
     * j - amount of columns.
     * 
     * Properties:
     * 1. (A^T)^T = A;
     * 2. (A + B)^T = A^T + B^T;
     * 3. (AB)^T = B^T * A^T;
     * 4. (n * A)^T = n * A^T wher n - any number;
     * 5. determinant A = determinant A^T.
     * @param {Number[][]} A Matrix A (m x n).
     * @returns {Number[][]} Transpose matrix A^T (n x m).
     */
    transpose: function(A) {
        var A_T = [];

        for (var i = 0; i < A[0].length; i++) {
            A_T[i] = [];
            for (var j = 0; j < A.length; j++) {
                A_T[j][i] = A[i][j]
            }
        }

        return A_T;
    },

    /**
     * Multiplies compatible matrices.
     * 
     * O(i * j * k) where:
     * 
     * i - amount of rows in A;
     * 
     * j - amount of columns in B;
     * 
     * k - amount of rows in B.
     * @param {Number[][]} A Matrix A.
     * @param {Number[][]} B Matrix B. 
     * @returns {Number[][]} Matrix AB.
     */
    multiply: function(A, B) {
        if (!A || !B || A.length === 0 || B.length === 0 || A[0].length !== B.length) {
            console.error(this.errors[0]);
            return false;
        }

        var AB = [];

        for (var i = 0; i < A.length; i++) {
            AB[i] = [];
            for (var j = 0; j < B[0].length; j++) {
                AB[i][j] = 0;
            }
        }
        
        for (var i = 0; i < A.length; i++) {
            for (var j = 0; j < B[0].length; j++) {
                var s = 0;
                for (var k = 0; k < B.length; k++) {
                    s += A[i][k] * B[k][j];
                }
                AB[i][j] = s;
            }
        }
        return AB;
    },

    /**
     * Solves the determinant of the source matrix by Gauss method with selection of the main element.
     * 
     * O(n^3)
     * @param {Number[][]} A Source matrix.
     * @param {Number} epsilon Error magnitude. Default: 1e-6.
     * @returns {Number} Determinant.
     */
    determinant: function(A, epsilon = 1e-6) {
        if (!A || A.length === 0 || A.length != A[0].length) {
            console.error(this.errors[1]);
            return false;
        }

        var determinant = 1;

        for (var i = 0; i < A.length; ++i) {
            var k = i;
            for (var j = i+1; j < A[0].length; ++j) {
                if (Math.abs(A[j][i]) > Math.abs(A[k][i])) {
                    k = j;
                }
            }

            if (Math.abs(A[k][i]) < epsilon) {
                determinant = 0;
                break;
            }

            [A[i], A[k]] = [A[k], A[i]];

            if (i != k) {
                determinant = -determinant;
            }

            determinant *= A[i][i];

            for (var j = i+1; j < A[0].length; ++j) {
                A[i][j] /= A[i][i];
            }

            for (var j = 0; j < A[0].length; ++j) {
                if (j != i && Math.abs(A[j][i]) > epsilon) {
                    for (var k = i+1; k < A[0].length; ++k) {
                        A[j][k] -= A[i][k] * A[j][i];
                    }
                }
            }
        }   
        return determinant;
    },

    /**
     * Computes the inverse matrix using Schultz's iterative method.
     * 
     * O(k * n^3) where:
     * 
     * k - amount of iterations;
     * 
     * n - size of the matrix A.
     * @param {Number[][]} A Source matrix.
     * @param {Number} epsilon Error magnitude. Default: 1e-6.
     * @param {Number} m Order of convergence of the method. Default 3.
     * @param {Number} max_iterations  Maximum number of iterations. Default: 1000.
     * @returns {Number[][]} Inversed matrix.
     */
    inverse: function(A, epsilon = 1e-6, m = 3, max_iterations = 1000) {
        if (!A || A.length === 0 || A[0].length == 0) {
            return false;
        }

        var E = this.identity(A.length);
        var A_norms = this.norma.first(A) * this.norma.infinity(A);
        var A_T = this.transpose(A);
        var U0 = [0];

        for (var i = 0; i < A.length; i++) {
            U0[i] = [];
            for (var j = 0; j < A.length; j++) {
                U0[i][j] = A_T[i][j] / A_norms;
            }
        }

        var k = 0;
        while (k <= max_iterations) {
            var psi = this.subtract(E, this.multiply(A, U0));
            var psi_norma = 0;
            for (var i = 0; i < psi.length; i++) {
                for (var j = 0; j < psi.length; j++) {
                    psi_norma += Math.abs(psi[i][j])**m;
                }
            }

            psi_norma = Math.sqrt(psi_norma);

            if (psi_norma <= epsilon) {
                return U0;
            } else {
                var psi_m_previous = psi;

                for (var i = 1; i <= m - 1; i++) {
                    psi_m = this.degree(psi_m_previous,i);
                    psi_m = this.add(psi_m_previous, psi_m);
                }
                var U1 = this.multiply(U0, this.add(E, psi_m_previous));
                U0 = U1;
                k++;
            }
        }
        console.error(this.errors[6]);
        return false;
    },
}
