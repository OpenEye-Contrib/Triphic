
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef NOTYET
int conn_table[25] = { 1 , 1 , 1 , 0 , 0 ,
		       1 , 1 , 0 , 1 , 1 ,
		       1 , 0 , 1 , 1 , 1 ,
		       0 , 1 , 1 , 1 , 1 ,
		       0 , 1 , 1 , 1 , 1 };
int conn_table[25] = { 1 , 1 , 1 , 1 , 1 ,
		       1 , 1 , 1 , 1 , 1 ,
		       1 , 1 , 1 , 1 , 1 ,
		       1 , 1 , 1 , 1 , 1 ,
		       1 , 1 , 1 , 1 , 1 };

int conn_table[25] = { 1 , 1 , 1 , 1 , 0 ,
		       1 , 1 , 0 , 1 , 1 ,
		       1 , 0 , 1 , 1 , 1 ,
		       1 , 1 , 1 , 1 , 0 ,
		       0 , 1 , 1 , 0 , 1 };

#endif

static const int MAX_CLIQUES = 100000;

#ifdef MALLOC
#include <dmalloc.h>
#endif

/************************************************************************/
/* this is the contents of file int_ind_shell_sortd.cc from DACLib,
   included here because this is a C file, and the C++ name mangler has
   confused the linker. There is no doubt a more elegant way of solving
   this problem. */
void int_ind_shell_sortd( int *pi_array , int *pi_index , int i_num ) {

  int   i , i_inc , i_v , i_value , j;

  i_inc = 1;
  do {
    i_inc *= 3;
    i_inc++;
  } while( i_inc <= i_num );

  pi_index--;
  do {
    i_inc /= 3;
    for( i = i_inc + 1 ; i <= i_num ; i++ ) {
      i_value = pi_array[pi_index[i]];
      i_v = pi_index[i];
      j = i;
      while( pi_array[pi_index[j - i_inc]] < i_value ) {
	pi_index[j] = pi_index[j - i_inc];
	j -= i_inc;
	if( j <= i_inc )
	  break;
      }
      pi_index[j] = i_v;
    }
  } while( i_inc > 1 );
}

/************************************************************************/
/* function to update the given set into the new set by transferring into
   it all nodes in the old set that are connected to the given node */
void update_sets( unsigned short node , int *old_sets , int old_not_end ,
		  int old_cand_end , int *new_sets , int *new_not_end ,
		  int *new_cand_end , char **adj_table ) {

  char     *conn_map = adj_table[node];

  int      i;

  *new_not_end = 0;
  for( i = old_not_end ; i ; i-- , old_sets++ )
    if( conn_map[*old_sets] ) {
      *new_sets = *old_sets;
      (*new_not_end)++;
      new_sets++;
    }
  *new_cand_end = *new_not_end;
  /* start from old_not_end + 1 to skip the current candidate which is
     at old_not_end */
  old_sets++;
  for( i = old_not_end + 1 ; i < old_cand_end ; i++ , old_sets++ )
    if( conn_map[*old_sets] ) {
      *new_sets = *old_sets;
      (*new_cand_end)++;
      new_sets++;
    }

}

/************************************************************************/
/* function select_candidate, to find the next node with which to extend the
   search. This is done by finding the first node in the cand_set which
   is not connected to the seed node. node_conns should contain the
   connection line for seed_node. The position of this node is returned */
int select_candidate( int seed_node , int *sets , int not_end , int cand_end ,
		      char **adj_table ) {

  char     *q;

  int      i , *p;

  /* find a node in cand_set not connected to seed_node, put it
     at the start of cand_set and return it */
  p = sets + not_end;
  q = adj_table[seed_node];
  for( i = not_end ; i < cand_end ; i++ , p++ )
    if( !q[*p] )
      break;

  return i;

}

/************************************************************************/
/* find the median clique size except that if the median equals the minimum,
   the median is incremented by 1, as the purpose is to reduce the number of
   cliques */
int find_median_clique_size( int num_cliques , int *clique_sizes ,
			     int max_cliques ) {

  int      i , cum_count , median_size , min_size , max_size , *size_counts;

  max_size = 0;
  min_size = INT_MAX;
  for( i = 0 ; i < num_cliques ; i++ ) {
    if( clique_sizes[i] > max_size )
      max_size = clique_sizes[i];
    if( clique_sizes[i] < min_size )
      min_size = clique_sizes[i];
  }

  size_counts = (int *) malloc( (size_t) ( max_size + 1 ) * sizeof( int ) );
  for( i = 0 ; i <= max_size ; i++ )
    size_counts[i] = 0;
  for( i = 0 ; i < num_cliques ; i++ )
    size_counts[clique_sizes[i]]++;

  cum_count = 0;
  for( i = 0 ; i <= max_size ; i++ ) {
    cum_count += size_counts[i];
    if( cum_count > num_cliques / 2 ) {
      median_size = i;
      break;
    }
  }

  if( median_size == min_size )
    median_size++;

#if DEBUG == 1
  printf( "Median clique size : %d (Max : %d and Min : %d)\n" ,
	  median_size , max_size , min_size );
#endif

  free( size_counts );
  return median_size;

}

/************************************************************************/
/* dump the all cliques that are smaller than the median clique size. The
   median clique size then becomes the minimum clique size, since we no
   longer want anything smaller than that. */
void dump_smallest_cliques( unsigned short **cliques , int *clique_sizes ,
			    int *num_cliques , int *max_cliques ,
			    int *min_clique_size ) {

  int      i , j , median_size;

  if( *max_cliques > MAX_CLIQUES )
    *max_cliques = MAX_CLIQUES;
  else
    *max_cliques = *num_cliques;

  median_size = find_median_clique_size( *num_cliques , clique_sizes ,
					 *max_cliques );

  if( median_size > *min_clique_size ) 
    *min_clique_size = median_size;

  /* ditch any cliques below the minimum size and shuffle the rest up */
  j = 0;
  for( i = 0 ; i < *num_cliques ; i++ ) {
    if( clique_sizes[i] < *min_clique_size ) {
      free( cliques[i] );
    } else {
      cliques[j] = cliques[i];
      clique_sizes[j] = clique_sizes[i];
      j++;
    }
  }
  *num_cliques = j;

#if DEBUG == 1
  printf( "New number of cliques : %d and new minimum clique size : %d\n" ,
	  *num_cliques , *min_clique_size );
#endif

}

/************************************************************************/
/* function to increase memory allocated for cliques. If requirements
   goes above MAX_CLIQUES or we run out of memory, dump all those with less
   than the median clique size rather than increasing memory */
void enlarge_clique_memory( unsigned short ***cliques , int **clique_sizes ,
			    int *num_cliques , int *max_cliques ,
			    int *min_clique_size ) {

  int      *pi_temp;

  unsigned short **ppi_temp;

  if( 0 == *max_cliques )
    *max_cliques = 10000;
  else
    (*max_cliques) *= 2;

  if( *max_cliques > MAX_CLIQUES )
    dump_smallest_cliques( *cliques , *clique_sizes , num_cliques ,
			   max_cliques , min_clique_size );

  /*  printf( "Trying to increase max_cliques to %d\n" , *max_cliques ); */
  ppi_temp = (unsigned short **) malloc( (size_t) (*max_cliques)
					 * sizeof( unsigned short * ) );
  if( NULL == ppi_temp ) {
    dump_smallest_cliques( *cliques , *clique_sizes , num_cliques ,
			   max_cliques , min_clique_size );
    ppi_temp = (unsigned short **) calloc( (size_t) (*max_cliques) ,
					   sizeof( unsigned short * ) );
    if( NULL == ppi_temp ) {
      /* it's really bad - so definitely give up */
      fprintf( stderr ,
	       "AWOOGA - out of memory in function enlarge_clique_memory\n\
of Bron-Kerbosch algorithm. You need to buy a bigger computer.\n" );
      fprintf( stderr , "Number of cliques : %d\n" , *num_cliques );
      exit( 1 );
    }
  }

  pi_temp = (int *) calloc( (size_t) (*max_cliques) , sizeof( int ) );
  if( NULL == pi_temp ) {
    fprintf( stderr ,
	     "AWOOGA - out of memory in function enlarge_clique_memory\n\
of Bron-Kerbosch algorithm. You need to buy a bigger computer.\n" );
    fprintf( stderr , "Number of cliques : %d\n" , *num_cliques );
    exit( 1 );
  }

  if( *num_cliques ) {
    memcpy( pi_temp , *clique_sizes , *num_cliques * sizeof( int ) );
    memcpy( ppi_temp , *cliques , *num_cliques * sizeof( int* ) );
  }

  if( *cliques ) free( *cliques );
  if( *clique_sizes ) free( *clique_sizes );
  *cliques = ppi_temp;
  *clique_sizes = pi_temp;

}

/************************************************************************/
/* find the node (min_node), which is the node in the not or candidate set
   with the least number of connections to other canditates. Also
   return the value of nd, the magic counter */
void find_min_connected_node( int *nd , int *min_node , int *min_node_pos ,
			      int *sets , int cand_end , int not_end ,
			      char **adj_table ) {

  char     *q;

  int      i , j , cand , pos , this_nd , *p , *r;

  *min_node = 0;
  cand = 0;
  *nd = cand_end;
  r = sets;
  for( i = 0 ; i < cand_end && *nd ; i++ , r++ ) {
    this_nd = 0;
    p = sets + not_end;
    q = adj_table[*r];
    for( j = not_end ; j < cand_end && this_nd < *nd ; j++ , p++ ) {
      if( !q[*p] ) {
	this_nd++;
	pos = j;
      }
    }

    if( this_nd < *nd ) {
      *min_node = *r;
      *nd = this_nd;
      if( i < not_end ) {
	*min_node_pos = pos;
      } else {
	*min_node_pos = i;
	cand = 1;
      }
    }
  }

  *nd += cand; /* increment magic count if selected node was a candidate */

}

/************************************************************************/
/* function extend - put the next possibility in the possible clique */
void extend( unsigned short *poss_clique , int depth , char **adj_table ,
	     int *in_sets , int in_cand_end , int in_not_end ,
	     unsigned short ***cliques , int **clique_sizes ,
	     int *num_cliques , int *max_cliques , int *min_clique_size ,
	     int largest_only ) {

  int      i , j , not_end , cand_end , min_node , min_node_pos , nd , *sets;

  /* if there aren't enough candidates left to make the clique as big as
     min_clique_size, bail out straightaway */
  if( depth + in_cand_end - in_not_end < *min_clique_size )
    return;

  /* get the memory for the local candidate and not lists, producing
     as little core store fragmentation as possible */
  sets = (int *) malloc( (size_t) ( in_cand_end + 1 ) * sizeof( int ) );
  if( NULL == sets ) {
    fprintf( stderr ,
	     "AWOOGA - out of memory in function extend of Bron-Kerbosch\n\
algorithm. You need to buy a bigger computer.\n" );
    fprintf( stderr , "Number of cliques : %d  depth %d\n" , *num_cliques ,
	     depth );
    exit( 1 );
  }

  /* find the fixed node, which is the node in the not or candidate set
     with the least number of connections to other candidates. Also
     return the value of nd, the magic counter */
  find_min_connected_node( &nd , &min_node , &min_node_pos , in_sets ,
			   in_cand_end , in_not_end , adj_table );

  /* this is the bron-kerbosch v2 criterion - only loop through nd times */

  for( i = nd ; i ; i-- ) {

    // do the swap with the new selected candidate
    j = in_sets[in_not_end];
    in_sets[in_not_end] = in_sets[min_node_pos];
    in_sets[min_node_pos] = j;

    poss_clique[depth] = in_sets[in_not_end];

    /* update the candidate and not sets by removing any in either not
       connected to the last member of poss_clique */
    update_sets( poss_clique[depth] , in_sets , in_not_end , in_cand_end ,
		 sets , &not_end , &cand_end , adj_table );

    /* if there are no candidates, need to do something */
    if( !cand_end ) {
      /* if there are also no nodes in the not set, it's a clique,
	 so store it if it's big enough */
      if( depth + 1 >= *min_clique_size ) {
	if( largest_only && depth + 1 > *min_clique_size ) {
	  /* need to ditch what we have so far as it's no longer big enough */
	  for( j = 0 ; j < *num_cliques ; j++ )
	    free( (*cliques)[j] );
	  *num_cliques = 0;
	  *min_clique_size = depth + 1;
	}
	(*clique_sizes)[*num_cliques] = depth + 1;
	(*cliques)[*num_cliques] =
	  (unsigned short *) malloc( (size_t) (depth+1) * sizeof( unsigned short ) );
	if( NULL == (*cliques)[*num_cliques] ) {
	  fprintf( stderr ,
		   "AWOOGA - out of memory in function extend of Bron-Kerbosch\n\
algorithm. You need to buy a bigger computer.\n" );
	  fprintf( stderr , "Number of cliques : %d  depth %d\n" , *num_cliques ,
		   depth + 1 );
	  exit( 1 );
	}
	
	memcpy( (*cliques)[*num_cliques] , poss_clique ,
		(depth + 1) * sizeof( unsigned short ) );
	(*num_cliques)++;
	if( *num_cliques == *max_cliques )
	  enlarge_clique_memory( cliques , clique_sizes , num_cliques ,
				 max_cliques , min_clique_size );
      }
    } else {
      if( not_end < cand_end )
	extend( poss_clique , depth + 1 , adj_table ,
		sets , cand_end , not_end , cliques , clique_sizes ,
		num_cliques , max_cliques , min_clique_size , largest_only );
    }
    /* put cand in not set and select a new one - this is achieved
       simply by incrementing in_not_end. */
    in_not_end++;

    /* find a candidate that is not connected to min_node and
       return it's position in in_sets */
    if( i )
      min_node_pos = select_candidate( min_node , in_sets , in_not_end ,
				       in_cand_end , adj_table );
  }

  free( sets );

}

/************************************************************************/
/* function bron_kerbosch to do the clique detection.

   Algorithm 457 - Finding all cliques of an undirected graph
   Communications of the ACM, vol 16, 575--577 (Sept 1973)

   Takes the correspondence graph in array conn_table. conn_table is a 2D
   array, giving for each row the nodes connected to it as numbers.  The
   1D array nconns gives the number of nodes connected to each row in the
   conn_table.
   Passes back a 2D array of the cliques in cliques, size of each clique in
   clique sizes and number of cliques in num_cliques.  All memory for these
   is allocated as required.
   NB - this implementation (which follows faithfully the recipe given in B
   & K's paper) requires that each node is flagged as connected to itself - it
   will most likely dump core at some point if this is not the case.

   This incarnation, mk4, works off an adjacency table rather than a connection
   table.  It should be faster, at least that's the hope.

   27th March 2002. Now, if min_clique_size comes in as -1, returns only the
   largest cliques i.e. all those of maximum size only.

*/
void bron_kerbosch( char **adj_table , int num_nodes ,
		    unsigned short ***cliques , int **clique_sizes ,
		    int *num_cliques , int *min_clique_size ) {

  int      i , cand_end , depth , largest_only , not_end , max_cliques , *sets;

  unsigned short  *poss_clique;

  max_cliques = 0;
  *num_cliques = 0;
  *clique_sizes = 0;
  *cliques = 0;

  for( i = 0 ; i < num_nodes ; i++ ) {
    if( !adj_table[i][i] ) {
      fprintf( stderr , "Major error in Bron-Kerbosch algorithm - requires\n");
      fprintf( stderr , "adjacency table diagonal elements to be non-zero.\n");
      exit( 1 );
    }
  }

  poss_clique =
    (unsigned short *) malloc( (size_t) num_nodes * sizeof( unsigned short ) );
  sets = (int *) malloc( (size_t) ( 2 * num_nodes ) * sizeof( int ) );

  enlarge_clique_memory( cliques , clique_sizes , num_cliques , &max_cliques ,
			 min_clique_size );

  /* fill the cand_set with all nodes */
  for( i = 0 ; i < num_nodes; i++ )
    sets[i] = i;
  cand_end = num_nodes;
  not_end = 0;

  depth = 0;

  if( -1 == *min_clique_size )
    largest_only = 1;
  else
    largest_only = 0;

  /* do the bron-kerbosch, which calls extend recursively till it's all over */
  extend( poss_clique , depth , adj_table , sets ,
	  cand_end , not_end , cliques , clique_sizes ,
	  num_cliques , &max_cliques , min_clique_size , largest_only );

  free( poss_clique );
  free( sets );

}
