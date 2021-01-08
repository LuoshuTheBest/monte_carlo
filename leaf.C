#include "leaf.H"
#include "pairtree.H"
#include "slist.H"

REAL CLeaf::m_distances[MAX_ROTAMER_SIZE][MAX_ROTAMER_SIZE];

// Switch the rotamer coordinates held by this leaf to different 
// precomputed values.
void CLeaf::changeRotamer(int rotIndex)
{
  assert(rotIndex < SIDECHAIN::m_aalist[getType()]->m_nRotamers);
  
  m_undoRotIndex = m_rotIndex;
  m_rotIndex = rotIndex;

  const ROTAMER & rot = 
    SIDECHAIN::m_aalist[getType()]->getRotamer(rotIndex, m_rotType);

  // Set the new coordinates as well as their BV and energy, which were
  // precomputed.
  m_positions = rot.m_positions;
  m_bv = rot.m_bv;
  m_energy = rot.m_energy;
}

// Compute all pairs of distances between the atoms of this leaf and
// the given leaf.
void CLeaf::computeDistances(CLeaf * pLeaf, const REAL rot[3][3],
			     const REAL trans[3])
{
  int size1 = SIDECHAIN::m_aalist[getType()]->m_size;
  int size2 = SIDECHAIN::m_aalist[pLeaf->getType()]->m_size;
  
  REAL cen[3], dist[3], vec1[3], vec2[3];

  // If this is a PRO backbone, we need to add the Cd atom because 
  // it is part of a group with the Ca atom for Electrostatic purposes
  if (getType() == BBP)
    {
      assert(getNext());
      assert(getNext()->getType() == PRO);
      MxVpV(vec1, m_rotate, getNext()->getPositions()[2], m_translate);
    }

  // If this is a PRO backbone, we need to add the Cd atom because 
  // it is part of a group with the Ca atom for Electrostatic purposes
  if (pLeaf->getType() == BBP)
    {
      assert(pLeaf->getNext());
      assert(pLeaf->getNext()->getType() == PRO);

      REAL temp[3];
      MxVpV(temp, pLeaf->m_rotate, pLeaf->getNext()->getPositions()[2], 
	    pLeaf->m_translate);
      MxVpV(vec2, rot, temp, trans); 
    }
  
  // Compute te distances between all pairs of atoms.
  for(int j = 0; j < size2; j++)
    {
      MxVpV(cen, rot, pLeaf->getPositions()[j], trans); 

      for (int i = 0; i < size1; i++)
	{	
	  VmV(dist, cen, getPositions()[i]);
	  m_distances[i][j] = Vlength2(dist);

	  // Add distances to the Cd of the second node (if type is BBP)
	  if (pLeaf->getType() == BBP)
	    {
	      VmV(dist, vec2, getPositions()[i]);
	      m_distances[i][size2] = Vlength2(dist);
	    }
	}

      // Add distances to the Cd of the first node (if type is BBP)
      if (getType() == BBP)
	{
	  VmV(dist, cen, vec1);
	  m_distances[size1][j] = Vlength2(dist);
	}
    }

  // Add distance between the Cd of the first and second nodes (both BBPs)
  if (getType() == BBP && pLeaf->getType() == BBP)
    {
      VmV(dist, vec2, vec1);
      m_distances[size1][size2] = Vlength2(dist);
    }
}

void CLeaf::computePairEnergy(CNode * pNode, const REAL rot[3][3],
			      const REAL trans[3], CTerm * term,
			      bool bSeparated)
{
  assert(pNode->isLeaf());
  CLeaf * pLeaf = (CLeaf*) pNode;

  // GLY node does not have any atoms. 
  // There is no interaction with it.
  if (getType() == GLY || pLeaf->getType() == GLY)
    return;

  // If the BVs are too far away, no need to do anything
  if (getBV()->computeDistance(pLeaf->getBV(), rot, trans) > CUTOFF_DISTANCE)
    {
      term->reset();
      return;
    }

  // Fill the pairwise distances matrix.
  computeDistances(pLeaf, rot, trans);

  REAL sum = 0.0;

  int diff = pLeaf->getIndex() - getIndex();

  // Compute all vdW terms.
  sum += CTerm::computeVdW(getType(), pLeaf->getType(),
  			    m_distances, diff);

  // Compute all elctrostatic terms
  sum += CTerm::computeElectrostatics(getType(), pLeaf->getType(),
  			      m_distances, diff);

  // Compute all Solvation terms.
  sum += CTerm::computeSolvation(getType(), pLeaf->getType(),
  				 m_distances, diff);

  // Store the energy sum at the corresponding leaf of the energytree.
  assert(term);
  term->set(sum);

  // Save a pointer to this leaf in case we need to undo the last move.
  m_undoPairs.push_back(term);
  
  return;
}

void CLeaf::computeSelfEnergy(CTerm * term)
{
  // Insert the precomputed energy of interaction between atoms inside
  // this leaf.
  term->set(m_energy);
  m_undoPairs.push_back(term);
  
  return;
}

bool CLeaf::findPairClash(CNode * pNode, const REAL rot[3][3],
			  const REAL trans[3], bool bSeparated)  
{
 assert(pNode->isLeaf());
 CLeaf * pLeaf = (CLeaf*) pNode;

 // GLY node does not have any atoms. 
 // There cannot be a clash with it.
 if (getType() == GLY || pLeaf->getType() == GLY)
   return false;
 
 int size1 = SIDECHAIN::m_aalist[getType()]->m_size;
 int size2 = SIDECHAIN::m_aalist[pLeaf->getType()]->m_size;
  
 REAL cen[3], dist[3];
 
 // Check all pairs of atoms for possible clash
 for(int j = 0; j < size2; j++)
   {
     MxVpV(cen, rot, pLeaf->getPositions()[j], trans); 
     int b = SIDECHAIN::m_aalist[pLeaf->getType()]->m_aTypes[j];

     for (int i = 0; i < size1; i++)
       {	
	 VmV(dist, cen, getPositions()[i]);
	 int a = SIDECHAIN::m_aalist[getType()]->m_aTypes[i];
	 
	 // Check exclusion list to see if this pair of atoms 
	 // should not be checked.
	 int ex = isExcluded(pLeaf->getIndex() - getIndex(), getType(), i, 
			     pLeaf->getType(), j);
	 if (ex != EXCLUDED &&
	     isStericClash(a, b, Vlength2(dist), ex == PAIR1_4))
	   {
	     return true;
	   }
       }
   }

 return false;
}

// Undo the changes to this leaf caused by the latest move.
void CLeaf::undo()
{
  if (m_nc & CNode::TRANSFORM)
    {
      m_angle = m_undoAngle;
      m_torsionE = m_undoTorsionE;
      McM(m_rotate, m_undoRotate);
    }

  if (m_nc & CNode::BOX)
    changeRotamer(m_undoRotIndex);
}
 
// Rotate the protein around the rotatable bond associated with 
// this leaf.
void CLeaf::rotate(REAL angle)
{ 
  // Store undo information.
  m_undoAngle = m_angle;
  McM(m_undoRotate, m_rotate);
  m_undoTorsionE = m_torsionE;

  // Change the angle and compute a new rotation matrix.
  m_angle += angle;
  compute_rotation(m_rotate, getJoint(), m_angle);

  // Recompute the torsion energy for this torsion angle.
  m_torsionE = compute_dihedral(getAngle(), getIndex() % 2);
}
