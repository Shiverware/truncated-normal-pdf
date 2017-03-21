
// http://picomath.org/javascript/phi.js.html
// CDF of (0,1;x)
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

    // Abramowitz and Stegun formula 7.1.26
    var t = 1.0 / (1.0 + p * x);
    var y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);

    return 0.5 * (1.0 + sign * y);
}

// PDF of (0,1;x)
function phi(x) {
    'use strict';
    var PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679,
        A = 1 / (Math.sqrt(2 * PI)),
        B = -1 / 2;

        return A * Math.exp( B * Math.pow(x, 2));
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

function truncatedMeanVar(a,b,mu,sigma) {
    'use strict';
    /* jshint newcap: false */
    var Alpha = (a - mu) / sigma;
    /* jshint newcap: false */
    var Beta = (b - mu) / sigma;
    /* jshint newcap: false */
    var PhiAlpha = Phi(Alpha);
    /* jshint newcap: false */
    var PhiBeta = Phi(Beta);
    var Z = PhiBeta - PhiAlpha;

    var tMuA = (phi(Alpha) - phi(Beta)) / Z;
    var tMu = mu + sigma * tMuA;

    var tVA = (Alpha * phi(Alpha) - Beta * phi(Beta)) / Z;
    var tVB = tMuA * tMuA;
    var tVariance = sigma * sigma * (1 + tVA - tVB);

    // console.log(a,b,mu,sigma, ' new ', tMu, tVariance);
    return [tMu, tVariance];
}

// console.log(tcdf([0.0, 1, 2.2, 5, 10, 15, 22], 0, 22, 5, 3));
// console.log(truncatedMeanVar(0, 22, 5, 3));
// console.log(tcdf([5.313409222051099, 5.313409222051099 - 2.7082698546589663, 5.313409222051099 + 2.7082698546589663], 0, 22, 5, 3));