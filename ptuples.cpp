/* Public Domain.  See the LICENSE file. */

/* Pythagorean Tuples generator                                */

/* A program to generate Pythagorean tuples.                   */
/* Uses a straightforward trial and error test if a given      */
/* tuple satisfies the equation.                               */
/* Pythagorean tuples are integer solutions to:                */
/*   a_1^2 + a_2^2 + ... = b^2                                 */

/* For N = 3, use the ptriples program which is much faster.   */

/* To compile, the GMP library needs to be already installed.  */
/* See https://gmplib.org                                      */
/* On linux, try:  gcc ptuples.cpp -lgmp -o ptuples            */

/* A great starting place for info is the Wikipedia page:      */
/* http://en.wikipedia.org/wiki/Pythagorean_triple             */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
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
void SaveToTuple( struct ttable*, long*, long, uint64_t*, unsigned long, mpz_t );
void MovePTuple( struct ttable*, mpz_t*, long, mpz_t );
void RemDupTuples( struct ttable* );
int CheckForDuplicateTuple( mpz_t*, mpz_t*, long );
long NumDupsAhead( struct ttable*, long );
int TupleIsPrimitive( mpz_t*, mpz_t, long );
void Cleanup_ttable( struct ttable* );
int ttable_tentry_cmpfunc( const void*, const void* );
int tentry_avalues_cmpfunc( const void*, const void* );

const unsigned int MAXB = 4294967294U;

int main( int argc, char * argv[] ) {

  if ( argc != 4 && argc != 5 ) {
    printf("\n");
    printf("For a_1^2 + a_2^2 + ... = b^2 :\n");
    printf("\n");
    printf("Usage: ptuples [-p] tuple_size b_min b_max\n\n\n");
    printf("Options:\n\n");
    printf("  -p -- primitive tuples only\n\n\n");
    printf("eg.  For all primitive Pythagorean quadruples from 100 to 500, try:\n\n");
    printf("ptuples -p 4 100 500\n\n");
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

  if ( mpz_cmp_ui( user_b_max, MAXB ) > 0 ) {
    printf("\nb_max must be <= %u.  Aborting.\n\n", MAXB );
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

void BuildNTuples( struct ttable* final_table, int DoOnlyPrimitives, mpz_t b_min, mpz_t b_max, long N ) {

  mpz_t b_min_sqr;
  mpz_init_set( b_min_sqr, b_min );
  mpz_mul( b_min_sqr, b_min_sqr, b_min_sqr );

  mpz_t b_max_sqr;
  mpz_init_set( b_max_sqr, b_max );
  mpz_mul( b_max_sqr, b_max_sqr, b_max_sqr );

  mpz_t Amax;
  mpz_init_set( Amax, b_max );
  mpz_sub_ui( Amax, Amax, 1 );

  unsigned long numsqrs = mpz_get_ui( Amax );
  uint64_t* sqrs = (uint64_t*) calloc( numsqrs, sizeof( uint64_t ) );
  if ( sqrs == NULL ) {
    printf("\nNot enough memory.  Exiting.\n");
    exit(-1);
  }

  long i;
  for ( i = 1; i <= numsqrs; i++ )
    sqrs[i-1] = i*i;

  const long numtumblers = N - 1;
  const long lasttumbler = numtumblers - 1;

  mpz_t* subtotaltumbler = (mpz_t*) calloc( numtumblers, sizeof( mpz_t ) );
  for ( i = 0; i < numtumblers; i++ )
    mpz_init( subtotaltumbler[i] );


  struct ttable tmp_table;
  tmp_table.count = 0;
  tmp_table.tuples = NULL;

  // initialize sqrindextumbler array
  long*  sqrindextumbler = (long*) calloc( numtumblers, sizeof( long ) );
  for ( i = 0; i < numtumblers; i++ )
    sqrindextumbler[i] = 0;

  int islessthan_b_min = 0;
  int isgreaterthan_b_max = 0;

  for ( i = 0; i >= 0; ) { // i is an index to the tumbler arrays.
    for ( ; i < lasttumbler; i++ ) {  // need to move i back to the last position
      mpz_set_ui( subtotaltumbler[i], sqrs[sqrindextumbler[i]] );
      if ( i > 0 )
        mpz_add( subtotaltumbler[i], subtotaltumbler[i], subtotaltumbler[i-1] );
    }

    // set last subtotal.  ie. The total of all squares of all a-values.
    mpz_add_ui( subtotaltumbler[i], subtotaltumbler[i-1], sqrs[sqrindextumbler[i]] );

    // optimization -- if first time, try to skip ahead
    if ( sqrindextumbler[i] == 0 ) {
      mpz_t tempZ;
      mpz_init( tempZ );

      mpz_sub( tempZ, b_min_sqr, subtotaltumbler[i-1] );

      if ( mpz_cmp_ui( tempZ, 1 ) >= 0 ) {
        mpz_sqrt( tempZ, tempZ );
        unsigned long sqrindex = mpz_get_ui( tempZ );
        if ( sqrindex > 0 )
          sqrindex--;
        if ( sqrindex >= numsqrs )
          sqrindex = numsqrs - 1;

        sqrindextumbler[i] = sqrindex;
        mpz_add_ui( subtotaltumbler[i], subtotaltumbler[i-1], sqrs[sqrindextumbler[i]] );
      }
      mpz_clear( tempZ );
    }

    islessthan_b_min = mpz_cmp( subtotaltumbler[i], b_min_sqr ) < 0;
    isgreaterthan_b_max = mpz_cmp( subtotaltumbler[i], b_max_sqr ) > 0;

    if ( !( islessthan_b_min || isgreaterthan_b_max ) )
      if ( mpz_perfect_square_p( subtotaltumbler[i] ) )
        SaveToTuple( &tmp_table, sqrindextumbler, numtumblers, sqrs, numsqrs, subtotaltumbler[i] );

    if ( isgreaterthan_b_max )  // already over the b_max_sqr limit -- skip to the end
      sqrindextumbler[i] = numsqrs;
    else
      sqrindextumbler[i]++;

    while ( i >= 0 && sqrindextumbler[i] >= numsqrs ) {
      sqrindextumbler[i] = 0;
      i--;
      if ( i >= 0 ) {
        if ( mpz_cmp( subtotaltumbler[i], b_max_sqr ) > 0 )  // already over the b_max_sqr limit -- skip to the end
         sqrindextumbler[i] = numsqrs;
        else
         sqrindextumbler[i]++;
      }
    }
  }

  RemDupTuples( &tmp_table );

  // move "tmp_table" tuples over to "final_table" table.
  for ( i = 0; i < tmp_table.count; i++ ) {
    if ( DoOnlyPrimitives )
      if ( !TupleIsPrimitive( tmp_table.tuples[i].a, tmp_table.tuples[i].b, tmp_table.tuples[i].a_count ) )
        continue;

    MovePTuple( final_table, tmp_table.tuples[i].a, tmp_table.tuples[i].a_count, tmp_table.tuples[i].b );
    tmp_table.tuples[i].a_count = 0;
    tmp_table.tuples[i].a = NULL;
  }

  Cleanup_ttable( &tmp_table );

  for ( i = 0; i < numtumblers; i++ )
    mpz_clear( subtotaltumbler[i] );
  free( subtotaltumbler );
  subtotaltumbler = NULL;

  free( sqrs );
  sqrs = NULL;

  mpz_clear( Amax );
  mpz_clear( b_max_sqr );
  mpz_clear( b_min_sqr );
}

void SaveToTuple( struct ttable* the_table, long* sqrindextumbler, long numtumblers, uint64_t* sqrs, unsigned long numsqrs, mpz_t bsqr ) {

  mpz_t* avalues = (mpz_t*)calloc( numtumblers, sizeof( mpz_t ) );

  long i;
  for ( i = 0; i <  numtumblers; i++ )
    mpz_init_set_ui( avalues[i], sqrindextumbler[i] + 1 );

  mpz_t b;
  mpz_init_set( b, bsqr );
  mpz_sqrt( b, b );

  MovePTuple( the_table, avalues, numtumblers, b );

  mpz_clear( b );
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

  long skip = 0;
  long i;
  for ( i = 0; i < the_table->count; i += skip + 1 ) {
    skip = NumDupsAhead( the_table, i ); 
    MovePTuple( &new_table, the_table->tuples[i].a, the_table->tuples[i].a_count, the_table->tuples[i].b );
    the_table->tuples[i].a_count = 0;
    the_table->tuples[i].a = NULL;
  }

  Cleanup_ttable( the_table );
  the_table->count = new_table.count;
  the_table->tuples = new_table.tuples;
}

long NumDupsAhead( struct ttable* the_table, long curindex ) {

  long DupCount = 0;

  long i = curindex;
  while ( i < the_table->count-1 &&
          CheckForDuplicateTuple( the_table->tuples[i].a, the_table->tuples[i+1].a, the_table->tuples[i].a_count ) ) {
    i++;
    DupCount++;
  }

  return DupCount;
}

int TupleIsPrimitive( mpz_t* avalues, mpz_t b, long acount ) {

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


