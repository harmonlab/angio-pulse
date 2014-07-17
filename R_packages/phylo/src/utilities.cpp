#include <Rcpp.h>
#include <vector>

/* 
 * GENERAL FUNCTIONS for dealing with an S3 'phylo' object from the R-package ape (Paradis E, J Claude, and K Strimmer 2004 [BIOINFORMATICS 20:289-290])
 * author: Jonathan M Eastman 01.11.2011
 */

/* INTERNAL C++ */
void sortedges(std::vector<double> const &unsortededges,
			   std::vector<double> &edges,
			   std::vector<int> const &des)
{
	std::vector<int>::size_type i,j;
	//int i, j;
	for(i=0; i<edges.size(); i++) {
		for(j=0; j<des.size(); j++){
			if(des.at(j)==(signed)i+1)
			{
				edges.at(i)=unsortededges.at(j);							 
			}
		}
	}
}

/* INTERNAL C++ */
void gatherdescendants(int const &node,
					   int const &root,
					   int const &endofclade,
					   std::vector<int> &TIPS,
					   std::vector<int> const &anc, 
					   std::vector<int> const &des,
					   int const &all)
{
	/* 
	 * node: current internal node
	 * root: most internal node
	 * endofclade: number of rows in the phylogenetic edge matrix (of the 'phylo' object)
	 * TIPS: vector in which to store all (tip) descendants of node (by node label)
	 * anc: first column of 'phylo' edge matrix (e.g., phy$edge[,1])
	 * des: second column of 'phylo' edge matrix (e.g., phy$edge[,2])
	 * all: whether to return all (all=1) or just tips (all=0)
	 */
	
	int i, eoc=endofclade, r=root, n=node, aa=all;

	for(i=0; i<eoc; i++) {
		if(anc.at(i) == n)
		{
			if(des.at(i)>r)
			{
				if(all==1)
				{
					TIPS.push_back(des.at(i));
				}
				
				gatherdescendants(des.at(i), r, eoc, TIPS, anc, des, aa);
				
			}
			else 
			{
				TIPS.push_back(des.at(i));
	
			}
		}
	}
}
	

/* 
 * FUNCTION TO return tip descendant given a node label and an S3 'phylo' object from the R-package ape (Paradis E, J Claude, and K Strimmer 2004 [BIOINFORMATICS 20:289-290])
 * author: Jonathan M Eastman 07.23.2011
 */

/* C++ | R INTERFACE; gather tip descendants main function */
RcppExport SEXP get_descendants (SEXP tree) 
{
	/* 
	 * tree: a list of elements 
	 * NODE: node for which descendants will be returned
	 * ROOT: most internal node
	 * ALL: whether to gather all (ALL=1) or just tips (ALL=0)
	 * ENDOFCLADE: rows in edge matrix
	 * ANC: first column of 'phylo' edge matrix (e.g., phy$edge[,1])
	 * DES: second column of 'phylo' edge matrix (e.g., phy$edge[,2])
	 */
	
	try {		
		/* call in parameters associated with 'phylo' object */
		Rcpp::List phylo(tree);
		
		int node = Rcpp::as<int>(phylo["NODE"]);
		int root = Rcpp::as<int>(phylo["ROOT"]);
		int all = Rcpp::as<int>(phylo["ALL"]);
		int endofclade =  Rcpp::as<int>(phylo["ENDOFCLADE"]);
		std::vector<int> anc=phylo["ANC"];
		std::vector<int> des=phylo["DES"];
		
		std::vector<int> TIPS;	
		if(all==0)
		{
			TIPS.reserve(root-1);
		}
		else 
		{
			TIPS.reserve(2*(root-1));	
		}
		
		/* get descendants */
		gatherdescendants(node, root, endofclade, TIPS, anc, des, all);

		/* PREPARE OUTPUT FOR R */
		return Rcpp::List::create(Rcpp::Named("TIPS",TIPS));
	
    } catch( std::exception &ex ) {		
		forward_exception_to_r( ex );
    } catch(...) { 
		::Rf_error( "C++ exception: unknown reason" ); 
    }
    return R_NilValue; 
}

/* 
 * FUNCTION TO GENERATE binary representation of all edges in tree from an S3 'phylo' object from the R-package ape (Paradis E, J Claude, and K Strimmer 2004 [BIOINFORMATICS 20:289-290])
 * author: Jonathan M Eastman 07.24.2011
 */

/* C++ | R INTERFACE; main function */
RcppExport SEXP binary_edges (SEXP tree) 
{
	/* 
	 * tree: a list of elements 
	 * NTIP: number of species in tree
	 * ROOT: most internal node
	 * ENDOFCLADE: rows in edge matrix
	 * MAXNODE: least internal internal node
	 * ANC: first column of 'phylo' edge matrix (e.g., phy$edge[,1])
	 * DES: second column of 'phylo' edge matrix (e.g., phy$edge[,2])
	 * EDGES: initialized binary matrix
	 */
	
	try {		
		/* call in parameters associated with 'phylo' object */
		Rcpp::List phylo(tree);
		
		int ntip = Rcpp::as<int>(phylo["NTIP"]);
		int root = Rcpp::as<int>(phylo["ROOT"]);
		int endofclade =  Rcpp::as<int>(phylo["ENDOFCLADE"]);
		int maxnode = Rcpp::as<int>(phylo["MAXNODE"]);
		std::vector<int> anc=phylo["ANC"];
		std::vector<int> des=phylo["DES"];
		std::vector<int> edges=phylo["EDGES"];

		int i,j;
		std::vector<int>::size_type k;
		int cur_node, cur_tip, cur_spp;
		
		std::vector<int> TIPS;	
		
		for(i=0; i<maxnode; i++)
			for(j=0; j<ntip; j++)
				edges[(ntip*i)+j]=0;
			
		
		/* collect tip descendants of each node */
		for(i=0; i<maxnode; i++){
			cur_node=i+1;
			if(cur_node<root)
			{
				cur_tip=cur_node;
				for(j=0; j<ntip; j++){
					cur_spp=j+1;
					if(cur_tip==cur_spp)
					{
						edges[(ntip*i)+j]=1;
					}
				}
			}
			else 
			{
				TIPS.reserve(root-1);
				gatherdescendants(cur_node, root, endofclade, TIPS, anc, des, 0);
				for(k=0; k<TIPS.size(); k++) {
					cur_tip=TIPS.at(k);
					for(j=0; j<ntip; j++){
						cur_spp=j+1;
						if(cur_tip==cur_spp)
						{
							edges[(ntip*i)+j]=1;
						}
					}
				}
				std::vector<int>().swap(TIPS);
			}
		}
		
		/* PREPARE OUTPUT FOR R */
		return Rcpp::List::create(Rcpp::Named("EDGES",edges));
		
    } catch( std::exception &ex ) {		
		forward_exception_to_r( ex );
    } catch(...) { 
		::Rf_error( "C++ exception: unknown reason" ); 
    }
    return R_NilValue; 
}
