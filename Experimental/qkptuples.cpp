/* Public Domain.  See the LICENSE file. */

/* QuicK Pythagorean Tuples generator                          */

/* A program to quickly generate SOME Pythagorean tuples.      */

/* Pythagorean tuples are integer solutions to:                */
/*   a_1^2 + a_2^2 + ... = b^2                                 */

/* To compile, the GMP library needs to be already installed.  */
/* See https://gmplib.org                                      */
/* On linux, try:  gcc qkptuples.cpp -lgmp -o qkptuples        */

/* A great starting place for info is the Wikipedia page:      */
/* http://en.wikipedia.org/wiki/Pythagorean_triple             */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>


struct tentry {
  long   a_count;
  mpz_t* a;
  mpz_t  b;
};

struct ttable {
long             count;
struct tentry*   tuples;
};


void BuildNTuples( struct ttable*, int, mpz_t, mpz_t, long );
void Build3Tuples( struct ttable*, int, mpz_t, mpz_t );
void AddPTriple( struct ttable*, mpz_t, mpz_t, mpz_t );
void MovePTuple( struct ttable*, mpz_t*, long, mpz_t );
void RemDupTuples( struct ttable* );
int CheckForDuplicateTuple( mpz_t*, mpz_t*, long );
int TupleIsPrimitive( mpz_t*, mpz_t, long );
long GetFirstBIndex( struct ttable*, mpz_t );
void Cleanup_ttable( struct ttable* );
int ttable_tentry_Bonly_cmpfunc( const void*, const void* );
int ttable_tentry_cmpfunc( const void*, const void* );
int tentry_avalues_cmpfunc( const void*, const void* );


int main( int argc, char * argv[] ) {

  if ( argc != 4 && argc != 5 ) {
    printf("\n");
    printf("For a_1^2 + a_2^2 + ... = b^2 :\n");
    printf("\n");
    printf("Usage: qkptuples [-p] tuple_size b_min b_max\n\n\n");
    printf("Options:\n\n");
    printf("  -p -- primitive tuples only\n\n\n");
    printf("eg.  For some primitive Pythagorean quadruples from 100 to 500, try:\n\n");
    printf("qkptuples -p 4 100 500\n\n");
    return 1;
  }

  int DoOnlyPrimitives = 0;
  if ( argc == 5 && strcmp( argv[1], "-p" ) == 0 )
    DoOnlyPrimitives = 1;

  long tuple_size = atol( argv[argc == 4 ? 1 : 2] );
  if ( tuple_size < 3 ) {
    printf("\ntuple size must be >= 3.  Aborting.\n\n");
    return 1;
  }

  mpz_t user_b_min;
  mpz_init_set_str( user_b_min,  argv[argc == 4 ? 2 : 3], 10 );

  mpz_t user_b_max;
  mpz_init_set_str( user_b_max,  argv[argc == 4 ? 3 : 4], 10 );

  if ( mpz_cmp_ui( user_b_min, 1 ) < 0 ) {
    printf("\nb_min must be >= 1.  Aborting.\n\n");
    mpz_clear( user_b_max );
    mpz_clear( user_b_min );
    return 1;
  }

  if ( mpz_cmp( user_b_min, user_b_max ) > 0 ) {
    printf("\nb_min must be <= b_max.  Aborting.\n\n");
    mpz_clear( user_b_max );
    mpz_clear( user_b_min );
    return 1;
  }

  struct ttable tuples;
  tuples.count = 0;
  tuples.tuples = NULL;

  BuildNTuples( &tuples, DoOnlyPrimitives, user_b_min, user_b_max, tuple_size );

  // print
  long i;
  for ( i = 0; i < tuples.count; i++ ) {
    printf("(");
    long j;
    for ( j = 0; j < tuples.tuples[i].a_count; j++ ) {
      if ( j > 0 )
        printf(",");
      gmp_printf("%Zd", tuples.tuples[i].a[j] );
    }
    gmp_printf(",%Zd)\n", tuples.tuples[i].b );
  }

  Cleanup_ttable( &tuples );

  mpz_clear( user_b_max );
  mpz_clear( user_b_min );

  return 0;
}

// Alorithm:
//   If N = 3, then just call Build3Tuples() and return.
//   Else
//    Build a 3 tuple table from b=5 to b=b_max
// A: Create a one larger tuple table by doing the following:
//    Check each a-value in the last built table to see if a
//    substitution can be made with an entry in the 3 tuple table.
//    Since a b^2 value in the 3 tuple table is equal to 2 squares,
//    we grow a tuple by one by substituting a single a-value with two
//    a-values from the 3-tuple table.  There may be multiple b-values
//    which match the a-value, so do it for all.
//    Repeat from Step A until we have the correct size table.
//    
//    NOTE: This will miss a lot of tuples, but should be a lot faster.
//
//    Some missed tuples:
//      2^2 + 2^2 + 1^2 = 3^2
//      1^2 + 1^2 + 1^2 + 1^2 = 2^2

void BuildNTuples( struct ttable* tuples, int DoOnlyPrimitives, mpz_t b_min, mpz_t b_max, long N ) {

  if ( N == 3 ) {
    Build3Tuples( tuples, DoOnlyPrimitives, b_min, b_max );
    qsort( tuples->tuples, tuples->count, sizeof(struct tentry), ttable_tentry_cmpfunc );
    return;
  }

  struct ttable threetuple;
  threetuple.count = 0;
  threetuple.tuples = NULL;

  mpz_t threetuple_min;
  mpz_init_set_ui( threetuple_min, 5 );

  Build3Tuples( &threetuple, 0, threetuple_min, b_max );
  qsort( threetuple.tuples, threetuple.count, sizeof(struct tentry), ttable_tentry_cmpfunc );

  struct ttable* onesmaller_table = NULL;
  struct ttable* beingbuilt_table = NULL;

  // Have to manually build the initial onesmaller table
  onesmaller_table = (struct ttable*) calloc( 1, sizeof( struct ttable ) );
  onesmaller_table->count = 0;
  onesmaller_table->tuples = NULL;

  long i;  
  for ( i = 0; i < threetuple.count; i++ ) {
    if ( mpz_cmp( threetuple.tuples[i].b, b_min ) >= 0 ) {
      mpz_t* avalues = (mpz_t*)calloc( 2, sizeof( mpz_t ) );
      mpz_init_set( avalues[0], threetuple.tuples[i].a[0] );
      mpz_init_set( avalues[1], threetuple.tuples[i].a[1] );
      MovePTuple( onesmaller_table, avalues, 2, threetuple.tuples[i].b );
    }
  }

  // iterate through all tuple sizes
  long tsize = 4;
  for ( tsize = 4; tsize <= N; tsize++ ) {
    if ( tsize > 4 ) {
      Cleanup_ttable( onesmaller_table );
      free( onesmaller_table );
      onesmaller_table = beingbuilt_table;
    }

    beingbuilt_table = (struct ttable*) calloc( 1, sizeof( struct ttable ) );
    beingbuilt_table->count = 0;
    beingbuilt_table->tuples = NULL;

    for ( i = 0; i < onesmaller_table->count; i++ ) {
      long j;
      for ( j = 0; j < onesmaller_table->tuples[i].a_count; j++ ) {
        long index = GetFirstBIndex( &threetuple, onesmaller_table->tuples[i].a[j] );
        if ( index >= 0 ) {
          do {
            long new_acount = tsize - 1;
            mpz_t* avalues = (mpz_t*)calloc( new_acount, sizeof( mpz_t ) );

            long k;
            for ( k = 0; k < new_acount; k++ )
              mpz_init( avalues[k] );

            for ( k = 0; k < j; k++ )
              mpz_set( avalues[k], onesmaller_table->tuples[i].a[k] );

            mpz_set( avalues[k], threetuple.tuples[index].a[0] );
            k++;

            mpz_set( avalues[k], threetuple.tuples[index].a[1] );
            k++;

            for ( ; k < new_acount; k++ )
              mpz_set( avalues[k], onesmaller_table->tuples[i].a[k-1] );

            MovePTuple( beingbuilt_table, avalues, new_acount, onesmaller_table->tuples[i].b );

            index++;
          } while ( index < threetuple.count && mpz_cmp( onesmaller_table->tuples[i].a[j], threetuple.tuples[index].b ) == 0 );
        }
      }
    }

    RemDupTuples( beingbuilt_table );
  }

  mpz_clear( threetuple_min );
  Cleanup_ttable( &threetuple );

// move "beingbuilt_table" tuples over to "tuples" table.
  for ( i = 0; i < beingbuilt_table->count; i++ ) {
    if ( DoOnlyPrimitives )
      if ( !TupleIsPrimitive( beingbuilt_table->tuples[i].a, beingbuilt_table->tuples[i].b, beingbuilt_table->tuples[i].a_count ) )
        continue;

    MovePTuple( tuples, beingbuilt_table->tuples[i].a, beingbuilt_table->tuples[i].a_count, beingbuilt_table->tuples[i].b );
    beingbuilt_table->tuples[i].a_count = 0;
    beingbuilt_table->tuples[i].a = NULL;
  }

  Cleanup_ttable( beingbuilt_table );

  qsort( tuples->tuples, tuples->count, sizeof(struct tentry), ttable_tentry_cmpfunc );
}

void Build3Tuples( struct ttable* tuples, int DoOnlyPrimitives, mpz_t b_min, mpz_t b_max ) {

  mpz_t working_b_min;
  if ( DoOnlyPrimitives )
    mpz_init_set( working_b_min, b_min );
  else
    mpz_init_set_ui( working_b_min, 1 );


  // We are going to use Euclid's formula:
  // For arbitrary positive integers m and n, and m > n,
  // a = m^2 - n^2, b = 2mn, c = m^2 + n^2

  // In addition, we will restrict GCD(m,n) = 1 and m-n be
  // an odd number to guarantee that the triple is primitive.

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

  // n can vary from 1 to no more than (b_max/2)^(1/2)
  mpz_t n_max;
  mpz_init_set( n_max, b_max );
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
    mpz_sub( m_min, working_b_min, n_squared );

    if ( mpz_cmp_ui( m_min, 1 ) < 0 )  // make sure mpz_sqrt() has a legal value
      mpz_set_ui( m_min, 1 );
    mpz_sqrt( m_min, m_min );
    mpz_sub_ui( m_min, m_min, 1 );  // subtract 1 just to be on the safe side

    // compute m_max
    mpz_sub( m_max, b_max, n_squared );
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
      if ( mpz_cmp( c, working_b_min ) < 0 )
        continue;
      if ( mpz_cmp( c, b_max ) > 0 )
        continue;

      if ( DoOnlyPrimitives )
        AddPTriple( tuples, a, b, c );
      else {
        // iterate through k in: (k*a)^2 + (k*b)^2 = (k*c)^2
        mpz_fdiv_q( k, b_min, c );

        for ( mpz_mul( kc, c, k ); mpz_cmp( kc, b_max ) <= 0; mpz_add_ui( k, k, 1 ), mpz_mul( kc, c, k ) ) {

          if ( mpz_cmp( kc, b_min ) < 0 )
            continue;

          mpz_mul( ka, a, k );
          mpz_mul( kb, b, k );

          AddPTriple( tuples, ka, kb, kc );
        }
      }
    }
  }

  mpz_clear( kc );
  mpz_clear( kb );
  mpz_clear( ka );
  mpz_clear( k );

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

  mpz_clear( working_b_min );
}


void AddPTriple( struct ttable* the_ttable, mpz_t a_0, mpz_t a_1, mpz_t b ) {

  mpz_t* avalues = (mpz_t*)calloc( 2, sizeof( mpz_t ) );
  mpz_init_set( avalues[0], a_0 );
  mpz_init_set( avalues[1], a_1 );

  MovePTuple( the_ttable, avalues, 2, b );
}


void MovePTuple( struct ttable* the_ttable, mpz_t* avalues, long a_count, mpz_t b ) {

  // allocate memory
  long index = 0;
  if ( the_ttable->count == 0 ) {
    the_ttable->count = 1;
    the_ttable->tuples = (struct tentry*) calloc( 1, sizeof(struct tentry) );
  }
  else {
    the_ttable->count++;
    the_ttable->tuples = (struct tentry*) realloc( the_ttable->tuples, sizeof(struct tentry) * the_ttable->count );
    index = the_ttable->count - 1;
    memset( &the_ttable->tuples[index], 0, sizeof(struct tentry) );
  }

  the_ttable->tuples[index].a = avalues;
  the_ttable->tuples[index].a_count = a_count;

  qsort( the_ttable->tuples[index].a, the_ttable->tuples[index].a_count, sizeof( mpz_t ), tentry_avalues_cmpfunc );

  mpz_init_set( the_ttable->tuples[index].b, b );
}

int CheckForDuplicateTuple( mpz_t* avalues1, mpz_t* avalues2, long count ) {

  int IsDuplicate = 1;

  long i;
  for ( i = 0; IsDuplicate && i < count; i++ )
    if ( mpz_cmp( avalues1[i], avalues2[i] ) != 0 )
      IsDuplicate = 0;

  return IsDuplicate;
}

void RemDupTuples( struct ttable* the_table ) {

  if ( the_table->count <= 1 )
    return;

  qsort( the_table->tuples, the_table->count, sizeof(struct tentry), ttable_tentry_cmpfunc );

  struct ttable new_table;
  new_table.count = 0;
  new_table.tuples = NULL;

  long i;
  for ( i = 0; i < the_table->count; i++ ) {
    int IsDuplicateTuple = 0;
    if ( i > 0 )
      IsDuplicateTuple = CheckForDuplicateTuple( the_table->tuples[i-1].a, the_table->tuples[i].a, the_table->tuples[i-1].a_count );

    if ( !IsDuplicateTuple ) {
      MovePTuple( &new_table, the_table->tuples[i].a, the_table->tuples[i].a_count, the_table->tuples[i].b );
      the_table->tuples[i].a_count = 0;
      the_table->tuples[i].a = NULL;
    }
  }

  Cleanup_ttable( the_table );
  the_table->count = new_table.count;
  the_table->tuples = new_table.tuples;
}

int TupleIsPrimitive( mpz_t* avalues, mpz_t b, long acount) {

  int IsPrimitive = 0;

  mpz_t gcd;
  mpz_init( gcd );

  if ( acount < 2 ) {
    printf("Internal Error:  TuplesIsPrimitive().  acount < 2.\n");
    return IsPrimitive;
  }

  mpz_gcd( gcd, avalues[0], avalues[1] );
  if ( mpz_cmp_ui( gcd, 1 ) == 0 )
    IsPrimitive = 1;

  long i;
  for ( i = 2; !IsPrimitive && i < acount; i++ ) {
    mpz_gcd( gcd, gcd, avalues[i] );
    if ( mpz_cmp_ui( gcd, 1 ) == 0 )
      IsPrimitive = 1;
  }

  if ( !IsPrimitive ) {
    mpz_gcd( gcd, gcd, b );
    if ( mpz_cmp_ui( gcd, 1 ) == 0 )
      IsPrimitive = 1;
  }

  mpz_clear( gcd );

  return IsPrimitive;
}

long GetFirstBIndex( struct ttable* the_table, mpz_t B ) {

  long index = -1;

  struct tentry entry;
  entry.a_count = 0;
  entry.a = NULL;
  mpz_init_set( entry.b, B );

  struct tentry* found_entry = (struct tentry*) bsearch( &entry, the_table->tuples, the_table->count,
                                sizeof( struct tentry ), ttable_tentry_Bonly_cmpfunc );

  mpz_clear( entry.b );

  if ( found_entry != NULL ) {
    index = found_entry - the_table->tuples;
    while ( index > 1 && mpz_cmp( the_table->tuples[index-1].b, B ) == 0 )
      index--;
  }

  return index;
}


// Free the memory allocated
void Cleanup_ttable( struct ttable* the_ttable ) {

  if ( the_ttable == NULL )
    return;

  long i;
  for ( i = 0; i < the_ttable->count; i++ ) {
    long j;
    for ( j = 0; j < the_ttable->tuples[i].a_count; j++ )
      mpz_clear( the_ttable->tuples[i].a[j] );
    free( the_ttable->tuples[i].a );
    the_ttable->tuples[i].a = NULL;
    mpz_clear( the_ttable->tuples[i].b );
  }

  if ( the_ttable->tuples != NULL ) {
    free( the_ttable->tuples );
    the_ttable->tuples = NULL;
  }

  the_ttable->count = 0;
}

int tentry_avalues_cmpfunc( const void* p1, const void* p2 ) {

  mpz_t*   avalue1 = (mpz_t*)p1;
  mpz_t*   avalue2 = (mpz_t*)p2;

  return mpz_cmp( *avalue1, *avalue2 );
}

int ttable_tentry_Bonly_cmpfunc( const void* p1, const void* p2 ) {

  struct tentry*   entry1 = (struct tentry*)p1;
  struct tentry*   entry2 = (struct tentry*)p2;

  return mpz_cmp( entry1->b, entry2->b );
}


int ttable_tentry_cmpfunc( const void* p1, const void* p2 ) {

  struct tentry*   entry1 = (struct tentry*)p1;
  struct tentry*   entry2 = (struct tentry*)p2;

  int cmpval = mpz_cmp( entry1->b, entry2->b );
  if ( cmpval == 0 ) {
    long i = 0;
    long a_count = entry1->a_count < entry2->a_count ? entry1->a_count : entry2->a_count;
    for ( i = 0; cmpval == 0 && i < a_count; i++ )
        cmpval = mpz_cmp( entry1->a[i], entry2->a[i] );
  }

  return cmpval;
}


