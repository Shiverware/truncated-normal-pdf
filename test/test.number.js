/* global describe, it, require */
'use strict';

// MODULES //

var // Expectation library:
	chai = require( 'chai' ),

	// Check whether an element is a finite number
	isFiniteNumber = require( 'validate.io-finite' ),

	// Check whether an element is `NaN`
	isnan = require( 'validate.io-nan' ),

	// Module to be tested:
	pdf = require( './../lib/number.js' );


// VARIABLES //

var expect = chai.expect,
	assert = chai.assert;


// TESTS //

describe( 'number pdf', function tests() {

	var validationData = require( './fixtures/number.json' ),
		data = validationData.data,
		expected = validationData.expected.map( function( d ) {
			return d === 'Inf' ? Infinity : d;
		}),
		a =  validationData.a,
		b = validationData.b,
		mu = validationData.mu,
		sigma = validationData.sigma;

	it( 'should export a function', function test() {
		expect( pdf ).to.be.a( 'function' );
	});

	it( 'should evaluate the truncated normal probability density function', function test() {
		var actual;
		for ( var i = 0; i < data.length; i++ ) {
			actual =  pdf( data[ i ], a, b, mu, sigma );
			if ( isFiniteNumber( actual ) && isFiniteNumber( expected[ i ] ) ) {
				assert.closeTo( actual, expected[ i ] , 1e-14 );
			}
		}
	});

	it( 'should return `0` if provided an `x` outside `[a,b]` as input', function test() {
		assert.strictEqual( pdf( a - 2, a, b, mu, sigma ), 0 );
		assert.strictEqual( pdf( b + 2, a, b, mu, sigma ), 0 );
	});

	it( 'should return `NaN` if provided `NaN` as input', function test() {
		assert.isTrue( isnan( pdf( NaN, a, b, mu, sigma ) ) );
	});

});
