/**
 * @file specialcoef.cpp
 * @author  C.Lehrenfeld
 * @date August 2015
 * @brief 
 */
#include "surffacetfes.hpp"
  
namespace bcip {
  using namespace ngcomp;

  template<ELEMENT_TYPE ET>
  ELEMENT_TYPE MyDummyFE<ET>::ElementType() const { return ET; }

  SurfaceFacetFESpace ::  
  SurfaceFacetFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    name="SurfaceFacetFESpace(l2surf)";
  }

  SurfaceFacetFESpace :: ~SurfaceFacetFESpace ()
  {
    ;
  }

  shared_ptr<FESpace> SurfaceFacetFESpace :: 
  Create (shared_ptr<MeshAccess> ma, const Flags & flags)
  {
    return make_shared<SurfaceFacetFESpace> (ma, flags, true);
  }

  void SurfaceFacetFESpace :: Update(LocalHeap & lh)
  {
    int nsel = ma->GetNSE();
    int nv = ma->GetNV();

    BitArray mark_vert(nv);
    mark_vert.Clear();

    vertdof.SetSize(nv);
    vertdof = -1;
    vertasso.SetSize(nv);
    vertdof = -1;
    ndof = 0;
    for (int elnr = 0; elnr < nsel; elnr++)
    {
      if (ma->GetSElIndex (elnr) == 1)
      {
	Ngs_Element ngel = ma->GetElement<1,BND> (elnr);
        for (int i = 0; i < 2; ++i)
        {
          if (!mark_vert.Test(ngel.vertices[i]))
            vertdof[ngel.vertices[i]] = ndof++;
          vertasso[ngel.vertices[i]] = elnr;
          mark_vert.Set(ngel.vertices[i]);
        }
      }
    }
    
  }

  const FiniteElement & SurfaceFacetFESpace :: GetSFE (int elnr, LocalHeap & lh) const
  {
    if (ma->GetSElIndex (elnr) == 1)
    {
      if (ma->GetDimension() == 2)
      {
	Ngs_Element ngel = ma->GetElement<1,BND> (elnr);
        return *new (lh) MyDummyFE<ET_SEGM>(vertasso[ngel.vertices[0]] == elnr,
                                            vertasso[ngel.vertices[1]] == elnr);
      }
      else
        throw Exception("no 3D yet");
    }
    else
      return *new (lh) DummyFE<ET_SEGM>();
  }
 
  const FiniteElement & SurfaceFacetFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    return *new (lh) DummyFE<ET_TRIG>();
    // Exception ("Volume elements not available for SurfaceFacetFESpace");
  }
 
  int SurfaceFacetFESpace :: GetNDof () const throw()
  {
    return ndof;
  }

  void SurfaceFacetFESpace :: GetSDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
    Ngs_Element ngel = ma->GetElement<1,BND> (elnr);
    for (int i = 0; i < 2; ++i)
      if (vertdof[ngel.vertices[i]] != -1)
        dnums.Append(vertdof[ngel.vertices[i]]);
  }
  
  void SurfaceFacetFESpace :: 
  GetDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize (0);
  }
  
  Table<int> * SurfaceFacetFESpace :: 
  CreateSmoothingBlocks ( int type) const
  {
    throw Exception ("no smoothing blocks yet");
    // int i, j, first;
    // Array<int> cnt(nel);
    // cnt = 0;
    // for (i = 0; i < nel; i++)
    //   cnt[i] = first_element_dof[i+1]-first_element_dof[i];
	
    // Table<int> & table = *new Table<int> (cnt);
    
    // for (i = 0; i < nel; i++)
    //   {
    //     first = first_element_dof[i];
    //     for (j = 0; j < cnt[i]; j++)
    //       table[i][j] = first+j;
    //   }
    // return &table;
  }


  void  SurfaceFacetFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  { dnums.SetSize(0); return; }
  
  void  SurfaceFacetFESpace ::GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  { dnums.SetSize(0); return; }
  
  void  SurfaceFacetFESpace ::GetFaceDofNrs (int fanr, Array<int> & dnums) const
  { dnums.SetSize(0); return; }
  
  void  SurfaceFacetFESpace ::GetInnerDofNrs (int elnr, Array<int> & dnums) const
  { dnums.SetSize(0); return; }
    




  // register FESpaces
  namespace surffacetfes_cpp
  {
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetFESpaceClasses().AddFESpace ("surffacet", SurfaceFacetFESpace::Create);
    }
    Init init;
  }
  
}//end of namespace bcip
