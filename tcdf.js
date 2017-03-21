
// http://picomath.org/javascript/phi.js.html
function Phi(x) {
    'use strict';
    // constants
    var a1 = 0.254829592;
    var a2 = -0.284496736;
    var a3 = 1.421413741;
    var a4 = -1.453152027;
    var a5 = 1.061405429;
    var p = 0.3275911;

    // Save the sign of x
    var sign = 1;
    if (x < 0) {
        sign = -1;
    }
    x = Math.abs(x) / Math.sqrt(2.0);

    // A&S formula 7.1.26
    var t = 1.0 / (1.0 + p * x);
    var y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);

    return 0.5 * (1.0 + sign * y);
}

// Return a function that can be used to calculate CDF value for truncated Normal
function partialTCDF(a, b, mu, sigma) {
    'use strict';
    /* jshint newcap: false */
    var PhiAlpha = Phi((a - mu) / sigma);
    /* jshint newcap: false */
    var PhiBeta = Phi((b - mu) / sigma);
    var Z = PhiBeta - PhiAlpha;
    var PhiEpsilonF = function (x) {
        return Phi((x - mu) / sigma);
    };

    return function (x) {
        if (x < a || x > b) {
            return 0;
        }
        /* jshint newcap: false */
        var PhiEpsilon = PhiEpsilonF(x);
        return (PhiEpsilon - PhiAlpha) / Z;
    };

}

// Process an array of values (x) for truncated normal with low (a), high(b)
function tcdf(x, a, b, mu, sigma) {
    'use strict';
    var len = x.length,
        fcn,
        i;
    var y = [];

    fcn = partialTCDF(a, b, mu, sigma);
    for (i = 0; i < len; i += 1) {
        if (typeof x[i] === 'number') {
            y[i] = fcn(x[i]);
        } else {
            y[i] = NaN;
        }
    }
    return y;
}

// console.log(tcdf([0.0, 1, 2.2, 5, 10, 15, 22], 0, 22, 5, 3));
