/* Public Domain.  See the LICENSE file. */

/* A simple program to generate Pythagorean triples.   */
/* ie.  Integer solutions to:  a^2 + b^2 = c^2         */

/* To compile, the GMP library needs to be already installed.    */
/* See https://gmplib.org                                        */
/* On linux, try:  gcc ptriples.cpp -lgmp -o ptriples            */

/* A great source of info is the Wikipedia page:      */
/* http://en.wikipedia.org/wiki/Pythagorean_triple    */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>


struct tentry {
  mpz_t a;
  mpz_t b;
  mpz_t c;
};

struct ttable {
long             count;
struct tentry*   triples;
};

void AddPTriple( struct ttable*, mpz_t, mpz_t, mpz_t );
void Cleanup_ttable( struct ttable* );
int ttable_entry_cmpfunc( const void*, const void* );


int main( int argc, char * argv[] ) {

  if ( argc != 3 && argc != 4 ) {
    printf("\n");
    printf("For a^2 + b^2 = c^2 :\n");
    printf("\n");
    printf("Usage: ptriples [-p] c_min c_max\n\n\n");
    printf("Options:\n\n");
    printf("  -p -- primitive triples only\n\n");
    return 1;
  }

  int DoOnlyPrimitives = 0;
  if ( argc == 4 && strcmp( argv[1], "-p" ) == 0 )
    DoOnlyPrimitives = 1;

  mpz_t user_c_min;
  mpz_init_set_str( user_c_min,  argv[argc == 3 ? 1 : 2], 10 );

  mpz_t user_c_max;
  mpz_init_set_str( user_c_max,  argv[argc == 3 ? 2 : 3], 10 );

  if ( mpz_cmp_ui( user_c_min, 1 ) < 0 ) {
    printf("\nc_min must be >= 1.  Aborting.\n\n");
    mpz_clear( user_c_max );
    mpz_clear( user_c_min );
    return 1;
  }

  if ( mpz_cmp( user_c_min, user_c_max ) > 0 ) {
    printf("\nc_min must be <= c_max.  Aborting.\n\n");
    mpz_clear( user_c_max );
    mpz_clear( user_c_min );
    return 1;
  }

  mpz_t working_c_min;
  if ( DoOnlyPrimitives )
    mpz_init_set( working_c_min, user_c_min );
  else
    mpz_init_set_ui( working_c_min, 1 );


  // We are going to use Euclid's formula:
  // For arbitrary positive integers m and n, and m > n,
  // a = m^2 - n^2, b = 2mn, c = m^2 + n^2

  // In addition, we will restrict GCD(m,n) = 1 and m-n be
  // an odd number to guarantee that the triple is primitive.

  struct ttable triples;
  triples.count = 0;
  triples.triples = NULL;

  mpz_t a;
  mpz_init( a );
  mpz_t b;
  mpz_init( b );
  mpz_t c;
  mpz_init( c );

  mpz_t n;
  mpz_init( n );
  mpz_t m;
  mpz_init( m );

  // n can vary from 1 to no more than (c_max/2)^(1/2)
  mpz_t n_max;
  mpz_init_set( n_max, user_c_max );
  mpz_cdiv_q_2exp( n_max, n_max, 1 );
  mpz_sqrt( n_max, n_max );

  mpz_t m_min;
  mpz_init( m_min );

  mpz_t m_max;
  mpz_init( m_max );

  mpz_t n_squared;
  mpz_init( n_squared );

  mpz_t m_squared;
  mpz_init( m_squared );

  mpz_t gcd;
  mpz_init( gcd );

  mpz_t tempZ;
  mpz_init( tempZ );

  mpz_t k;
  mpz_init( k );

  mpz_t ka;
  mpz_init( ka );
  mpz_t kb;
  mpz_init( kb );
  mpz_t kc;
  mpz_init( kc );

  // iterate through n
  for ( mpz_set_ui( n, 1 );  mpz_cmp( n, n_max ) <= 0; mpz_add_ui( n, n, 1 ) ) {
    mpz_mul( n_squared, n, n );

    // compute m_min
    mpz_sub( m_min, working_c_min, n_squared );

    if ( mpz_cmp_ui( m_min, 1 ) < 0 )  // make sure mpz_sqrt() has a legal value
      mpz_set_ui( m_min, 1 );
    mpz_sqrt( m_min, m_min );
    mpz_sub_ui( m_min, m_min, 1 );  // subtract 1 just to be on the safe side

    // compute m_max
    mpz_sub( m_max, user_c_max, n_squared );
    if ( mpz_cmp_ui( m_max, 1 ) < 0 )  // make sure mpz_sqrt() has a legal value
      mpz_set_ui( m_max, 1 );
    mpz_sqrt( m_max, m_max );

    // calc first value of m
    if ( mpz_cmp( n, m_min ) < 0 ) {
      mpz_set( m, m_min );
      mpz_sub( tempZ, m, n );
      if ( mpz_divisible_ui_p( tempZ, 2 ) )
        mpz_add_ui( m, m, 1 );
    }
    else {
      mpz_set( m, n );
      mpz_add_ui( m, m, 1 );
    }

    // iterate through m
    for ( ; mpz_cmp( m, m_max ) <= 0; mpz_add_ui( m, m, 2 ) ) {

      // generate a primitive (a,b,c)
      mpz_gcd( gcd, m, n );
      if ( mpz_cmp_ui( gcd, 1 ) != 0 )
        continue;

      mpz_mul( m_squared, m, m );

      mpz_sub( a, m_squared, n_squared );

      mpz_mul( b, m, n );
      mpz_mul_ui( b, b, 2 );

      mpz_add( c, m_squared, n_squared );

      // check if primitive is outside our working range
      if ( mpz_cmp( c, working_c_min ) < 0 )
        continue;
      if ( mpz_cmp( c, user_c_max ) > 0 )
        continue;

      if ( DoOnlyPrimitives )
        AddPTriple( &triples, a, b, c );
      else {
        // iterate through k in: (k*a)^2 + (k*b)^2 = (k*c)^2
        mpz_fdiv_q( k, user_c_min, c );

        for ( mpz_mul( kc, c, k ); mpz_cmp( kc, user_c_max ) <= 0; mpz_add_ui( k, k, 1 ), mpz_mul( kc, c, k ) ) {

          if ( mpz_cmp( kc, user_c_min ) < 0 )
            continue;

          mpz_mul( ka, a, k );
          mpz_mul( kb, b, k );

          AddPTriple( &triples, ka, kb, kc );
        }
      }
    }
  }

  mpz_clear( kc );
  mpz_clear( kb );
  mpz_clear( ka );
  mpz_clear( k );

  qsort( triples.triples, triples.count, sizeof(struct tentry), ttable_entry_cmpfunc );

  // print
  long i;
  for ( i = 0; i < triples.count; i++ )
    gmp_printf("(%Zd,%Zd,%Zd)\n", triples.triples[i].a, triples.triples[i].b, triples.triples[i].c );

  mpz_clear( tempZ );
  mpz_clear( gcd );
  mpz_clear( m_squared );
  mpz_clear( n_squared );
  mpz_clear( m_max );
  mpz_clear( m_min );
  mpz_clear( n_max );
  mpz_clear( m );
  mpz_clear( n );
  mpz_clear( c );
  mpz_clear( b );
  mpz_clear( a );

  Cleanup_ttable( &triples );

  mpz_clear( working_c_min );
  mpz_clear( user_c_max );
  mpz_clear( user_c_min );

  return 0;
}

// Add an entry to TABLE2
void AddPTriple( struct ttable* the_ttable, mpz_t a, mpz_t b, mpz_t c ) {

  // allocate memory
  long index = 0;
  if ( the_ttable->count == 0 ) {
    the_ttable->count = 1;
    the_ttable->triples = (struct tentry*) calloc( 1, sizeof(struct tentry) );
  }
  else {
    the_ttable->count++;
    the_ttable->triples = (struct tentry*) realloc( the_ttable->triples, sizeof(struct tentry) * the_ttable->count );
    index = the_ttable->count - 1;
    memset( &the_ttable->triples[index], 0, sizeof(struct tentry) );
  }

  mpz_init_set( the_ttable->triples[index].a, mpz_cmp( a, b ) < 0 ? a : b );
  mpz_init_set( the_ttable->triples[index].b, mpz_cmp( a, b ) < 0 ? b : a );
  mpz_init_set( the_ttable->triples[index].c, c );
}

// Free the memory allocated
void Cleanup_ttable( struct ttable* the_ttable ) {

  if ( the_ttable == NULL )
    return;

  long i;
  for ( i = 0; i < the_ttable->count; i++ ) {
    mpz_clear( the_ttable->triples[i].a );
    mpz_clear( the_ttable->triples[i].b );
    mpz_clear( the_ttable->triples[i].c );
  }

  if ( the_ttable->triples != NULL ) {
    free( the_ttable->triples );
    the_ttable->triples = NULL;
  }

  the_ttable->count = 0;
}

int ttable_entry_cmpfunc( const void* p1, const void* p2 ) {

  struct tentry*   entry1 = (struct tentry*)p1;
  struct tentry*   entry2 = (struct tentry*)p2;

  int cmpval = mpz_cmp( entry1->c, entry2->c );
  if ( cmpval == 0 )
      cmpval = mpz_cmp( entry1->b, entry2->b );
  if ( cmpval == 0 )
      cmpval = mpz_cmp( entry1->c, entry2->c );

  return cmpval;
}


