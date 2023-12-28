const scale_factor = 2;
var scale = 1;
var canvas_scale = 1;
var canvas = document.getElementById('canvas');
var context = canvas.getContext('2d');
const canvas_width = document.documentElement.scrollWidth;
const canvas_height = document.documentElement.scrollHeight;
const canvas_inner_width = canvas_width * scale_factor;
const canvas_inner_height = canvas_height * scale_factor;
const buffer = context.getImageData(0, 0, canvas.width, canvas.height);
const step = buffer.width * 4;
const dpi = 1/2.54;
const ratio = context.miterLimit*dpi*scale_factor;
canvas.width = canvas_inner_width;
canvas.height = canvas_inner_height;
canvas.style.width = canvas_width + 'px';
canvas.style.height = canvas_height + 'px';

const rad_deg = 180/Math.PI;
const deg_rad = Math.PI/180;

var position = {
    x: 0, 
    y: 0
};

var current_area = {
    x0: 0, 
    y0: 0,
    x1: canvas_width,
    y1: canvas_height
};

/**
 * @class
 * @param {Number} x Coordinate x.
 * @param {Number} y Coordinate y.
 */
class NODE {
    constructor(x, y) {
        this.x = x;
        this.y = y;
    }
}

/**
 * @class
 * @param {NODE} node_0 Start point of the line.
 * @param {NODE} node_1 End point of the line.
 */
class LINE {
    constructor(node_0, node_1) {
        this.x0 = node_0.x;
        this.y0 = node_0.y;
        this.x1 = node_1.x;
        this.y1 = node_1.y;
    }
}

class ELEMENT {
    constructor(node_1, node_2, fixing, E, A, n) {
        const line = new LINE(node_1, node_2);
        this.x0 = line.x0;
        this.y0 = line.y0;
        this.x1 = line.x1;
        this.y1 = line.y1;
        this.fixing = fixing;
        this.E = E;
        this.A = A;
        this.n = n;
    }

    length() {
        return Math.sqrt((this.x1 - this.x0) ** 2 + (this.y1 - this.y0) ** 2);
    }

    angle() {
        return Math.atan((this.y1 - this.y0) / (this.x1 - this.x0)) * rad_deg;
    }

    cos() {
        return (this.x1 - this.x0) / this.length();
    }

    sin() {
        return (this.y1 - this.y0) / this.length();
    }
}

var DRAW = {
    /**
     * Takes into account the constancy of line thickness during scaling
     * @param {Boolean} zooming True - the line is scaled, false - the line has a constant size.
     * @returns {Number} Scale factor.
     */
    multiplier: function(zooming) {
        if (zooming === true) {
            return ratio;
        } else {
            return ratio/scale;
        }
    },

    /**
     * Represents lines in accordance with GOST 2.303-68 (government standard).
     */
    line: {
        /**
         * @param {LINE} line An instance of the LINE class representing the start and end points of a line.
         * @param {Number[]} color Default - black.
         * @param {Number} width Default - 1mm.
         * @param {Boolean} zooming - Scale factor, default true.
         */
        main: function(line, color=[0,0,0,255], width = 1, zooming = true) {
            context.beginPath();
            context.lineWidth = width*DRAW.multiplier(zooming);
            context.lineJoin = 'round';
            context.lineCap = 'round';
            context.strokeStyle = `rgba(${color[0]},${color[1]},${color[2]},${color[3]})`;

            context.moveTo((line.x0*ratio - current_area.x0), (line.y0*ratio - current_area.y0));
            context.lineTo((line.x1*ratio - current_area.x0), (line.y1*ratio - current_area.y0));
            context.stroke();
            context.closePath();
        },

        /**
         * @param {LINE} line An instance of the LINE class representing the start and end points of a line.
         * @param {Number[]} color Default - black.
         * @param {Number} width Default - 0,5mm.
         * @param {Boolean} zooming - Scale factor, default true.
         */
        thin: function(line, color=[0,0,0,255], width = 0.5, zooming = true) {
            context.beginPath();
            context.lineWidth = width*DRAW.multiplier(zooming);
            context.lineJoin = 'round';
            context.lineCap = 'round';
            context.strokeStyle = `rgba(${color[0]},${color[1]},${color[2]},${color[3]})`;
            context.moveTo((canvas_inner_width/2 + line.x0*ratio), (canvas_inner_height/2 - line.y0*ratio));
            context.lineTo((canvas_inner_width/2 + line.x1*ratio), (canvas_inner_height/2 - line.y1*ratio));
            context.stroke();
            context.closePath();
        },

        /**
         * @param {LINE} line An instance of the LINE class representing the start and end points of a line.
         * @param {Number[]} color Default - black.
         * @param {Number} width Default - 0,5mm.
         * @param {Boolean} zooming - Scale factor, default true.
         */
        dashed: function(line, color=[0,0,0,255], width = 0.5, zooming = true) {
            context.beginPath();
            context.lineWidth = width*DRAW.multiplier(zooming);
            context.lineJoin = 'round';
            context.lineCap = 'round';
            context.strokeStyle = `rgba(${color[0]},${color[1]},${color[2]},${color[3]})`;
            context.moveTo((canvas_inner_width/2 + line.x0*ratio), (canvas_inner_height/2 - line.y0*ratio));
            context.lineTo((canvas_inner_width/2 + line.x1*ratio), (canvas_inner_height/2 - line.y1*ratio));
            context.setLineDash([2*ratio, 2*ratio]);
            context.stroke();
            context.closePath();
        },
    },

    /**
     * @param {LINE} line An instance of the LINE class representing the start and end points of a line.
     * @param {Number[]} color Default - black.
     * @param {Number} width Default - 0,5mm.
     * @param {Number} arrow_length The default is 2.5mm according to the standard. 
     * @param {Boolean} zooming - Scale factor, default false.
     */
    arrow: function(line, color = [0, 0, 0, 255], width = 0.25, arrow_length = 5, zooming = false) {
        var x0 = line.x0;
        var y0 = line.y0;
        var x1 = line.x1;
        var y1 = line.y1;
        var length = ((x1-x0)**2+(y1-y0)**2)**(1/2);
        var tan_alpha = (y1 - y0) / (x1 - x0);
        var cos_alpha = (x1-x0) / length;
        var sin_alpha = (y1-y0) / length;
        var tan_beta = -1 / tan_alpha;
        var beta = Math.atan(tan_beta) * rad_deg;
        var a = arrow_length * Math.tan(10 * deg_rad);
        var base = [
            (length-arrow_length)*cos_alpha + x0,
            (length-arrow_length)*sin_alpha + y0
        ];
        var points = [
            base[0] + a * Math.cos(beta * deg_rad),
            base[1] + a * Math.sin(beta * deg_rad),
            base[0] - a * Math.cos(beta * deg_rad),
            base[1] - a * Math.sin(beta * deg_rad)
        ];

        context.beginPath();
        context.lineWidth = width*DRAW.multiplier(zooming);
        context.lineJoin = 'round';
        context.lineCap = 'round';
        context.strokeStyle = `rgba(${color[0]},${color[1]},${color[2]},${color[3]})`;
        context.moveTo((canvas_inner_width / 2 + x0 * ratio), (canvas_inner_height / 2 - y0 * ratio));
        context.lineTo((canvas_inner_width / 2 + x1 * ratio), (canvas_inner_height / 2 - y1 * ratio));
        context.stroke();
        context.closePath();

        context.beginPath();
        context.lineWidth = width*DRAW.multiplier(zooming);
        context.lineJoin = 'miter';
        context.lineCap = 'round';
        context.strokeStyle = `rgba(${color[0]},${color[1]},${color[2]},${color[3]})`;
        context.fillStyle = `rgba(${color[0]},${color[1]},${color[2]},${color[3]})`;
        context.moveTo((canvas_inner_width / 2 + points[0] * ratio), (canvas_inner_height / 2 - points[1] * ratio));
        context.lineTo((canvas_inner_width / 2 + x1 * ratio), (canvas_inner_height / 2 - y1 * ratio));
        context.lineTo((canvas_inner_width / 2 + points[2] * ratio), (canvas_inner_height / 2 - points[3] * ratio));
        context.closePath();
        context.stroke();
        context.fill();
    },

    /**
     * 
     * @param {NODE} node An instance of the NODE class representing coordinates of the piont.
     * @param {Number} size Length of the arrows.
     */
    coordinate_system: function(node, size) {
        DRAW.arrow(new LINE(new NODE(node.x,node.y), new NODE(node.x+size,node.y)), [0,255,0,255], 0.25, 2.5, false)
        DRAW.arrow(new LINE(new NODE(node.x,node.y), new NODE(node.x,node.y+size)), [255,0,0,255], 0.25, 2.5, false)
    },

    circle: {
        /** 
         * @param {NODE} node An instance of the NODE class representing coordinates of the piont.
         * @param {Number} diameter Size of the diameter in mm.
         * @param {Number[]} color Default - black.
         * @param {Number} width Default 0.25mm.
         * @param {Boolean} zooming - Scale factor, default true.
         */
        fill: function(node, diameter, color = [0, 0, 0, 255], width = 0.25, zooming = true) {
            context.beginPath();
            context.lineWidth = width*DRAW.multiplier(zooming);
            context.fillStyle = `rgba(${color[0]},${color[1]},${color[2]},${color[3]})`;
            context.arc((canvas_inner_width/2 + node.x*ratio), (canvas_inner_height/2 - node.y*ratio), (diameter/2)*DRAW.multiplier(zooming), 0, 2*Math.PI);
            context.fill();
            context.closePath();
        },

        /** 
         * @param {NODE} node An instance of the NODE class representing coordinates of the piont.
         * @param {Number} diameter Size of the diameter in mm.
         * @param {Number[]} color Default - black.
         * @param {Number} width Default 0.25mm.
         * @param {Boolean} zooming - Scale factor, default true.
         */
        stroke: function(node, diameter, color = [0, 0, 0, 255], width = 0.25, zooming = true) {
            context.beginPath();
            context.lineWidth = width*DRAW.multiplier(zooming);
            context.strokeStyle = `rgba(${color[0]},${color[1]},${color[2]},${color[3]})`;
            context.arc((canvas_inner_width/2 + node.x*ratio), (canvas_inner_height/2 - node.y*ratio), (diameter/2)*DRAW.multiplier(zooming), 0, 2*Math.PI);
            context.stroke();
            context.closePath();
        },
    },
}

var CANVAS = {
    update: function() {
        context.putImageData(buffer, 0, 0);
    },

    /**
     * Represents all elelements in the canvas.
     */
    draw: function() {
        if (context === null) {return;}

        DRAW.coordinate_system(new NODE(0,0), 20)
        DRAW.circle.fill(new NODE(0,0),1,[0,0,0,255],false);
        DRAW.circle.stroke(new NODE(0,0),2,[0,0,0,255]);
    },

    /**
     * Scales the canvas when the mouse wheel is scrolled relative to the current cursor position.
     */
    scale: function() {
        canvas.addEventListener("wheel", onmousewheel, false);
        canvas.addEventListener("DOMMouseScroll", onmousewheel, false);
    
        function onmousewheel(event) {
            var mouse_x = event.clientX*2;
            var mouse_y = event.clientY*2;
            var delta = event.type === "wheel" ? event.wheelDelta : -event.detail;
            var increment = delta > 0 ? 1.1 : 0.9;
            scale *= increment;

            position.x = mouse_x;
            position.y = mouse_y;

            current_area.x0 -= mouse_x;
            current_area.y0 -= mouse_y;
            current_area.x1 -= mouse_x;
            current_area.y1 -= mouse_y;
            current_area.x0 *= increment;
            current_area.y0 *= increment;
            current_area.x1 *= increment;
            current_area.y1 *= increment;
            current_area.x0 += mouse_x;
            current_area.y0 += mouse_y;
            current_area.x1 += mouse_x;
            current_area.y1 += mouse_y;

            function redraw() {
                canvas_scale = -(current_area.x0 - current_area.x1) / canvas_width;
                context.setTransform(1, 0, 0, 1, 0, 0);
                context.clearRect(0, 0, canvas_inner_width, canvas_inner_height);
                context.setTransform(canvas_scale, 0, 0, canvas_scale, current_area.x0, current_area.y0);
                CANVAS.draw();
                requestAnimationFrame(() => redraw());
            }
            redraw();
            event.preventDefault();
        }
    },
}
CANVAS.draw();
CANVAS.scale();

document.addEventListener("mousemove", event => {
    position.x = event.clientX*2;
    position.y = event.clientY*2;
})

var MATRIX = {
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
            return MATRIX.selection_sort(first_norma)[0][0];
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
            return MATRIX.selection_sort(infinity_norma)[0][0];
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
                A_T[i][j] = A[j][i];
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
