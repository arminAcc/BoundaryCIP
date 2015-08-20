/**
 * @file specialcoef.hpp
 * @author  C.Lehrenfeld
 * @date August 2015
 * @brief
 */

#pragma once 

#include <comp.hpp>



namespace bcip {


  using namespace ngcomp;



  /**
     a placeholder finite element - 
   */
  template <ELEMENT_TYPE ET>
  class MyDummyFE : public FiniteElement
  {
  public:
    bool pos[2] = {false,false};
    /* INLINE */ MyDummyFE (bool leftpos, bool rightpos) : FiniteElement(2, 0)
    {
      pos[0] = leftpos;
      pos[1] = rightpos;
    }
    HD virtual ELEMENT_TYPE ElementType() const;
  };


  
  class NGS_DLL_HEADER SurfaceFacetFESpace : public FESpace
  {
  protected:
  
    // Level
    int level;
    Array<int> vertdof;
    Array<int> vertasso;
    int ndof;


    // if order is relative to mesh order 
    bool var_order = false;
    // variable order is set to mesh_order + rel_order 
    int rel_order = 0;

  public:

    SurfaceFacetFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false);
    ///
    virtual ~SurfaceFacetFESpace ();

    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags);

    virtual string GetClassName () const
    {
      return "SurfaceFacetFESpace";
    }

    ///
    virtual void Update(LocalHeap & lh);
    /// 
    //virtual void UpdateDofTables();
    ///
    virtual int GetNDof () const throw();
    ///
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    ///
    virtual const FiniteElement & GetSFE (int elnr, LocalHeap & lh) const;
    ///
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;
  
    virtual Table<int> * CreateSmoothingBlocks ( int type = 0) const;

    /// 
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const;
    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const;
    virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const;
    virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const;

    virtual bool VarOrder() const { return var_order; } 
    virtual int GetRelOrder() const { return rel_order; }   

  };

}//end of namespace bcip



