/*
 *  $RCSfile$
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright (C) 1997-2005  The Chemistry Development Kit (CDK) project
 *
 *  Contact: cdk-devel@lists.sourceforge.net
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2.1
 *  of the License, or (at your option) any later version.
 *  All we ask is that proper credit is given for our work, which includes
 *  - but is not limited to - adding the above copyright notice to the beginning
 *  of your source code files, and to any copyright notice that you may distribute
 *  with programs based on this work.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */
package org.openscience.cdk.graph;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.graph.AtomContainerPermutor;;


/**
 * This class allows to iterate trough the set of all possible
 * permutations of the bond order in a given atom container. 
 * It allows to check for a dependency of algorithm results 
 * on the bond order.
 *
 * The permutation code here is based on a pseudo code example 
 * on a tutorial site created and maintained by Phillip P. Fuchs:
 * http://www.geocities.com/permute_it/pseudo2.html
 * 
 *@author         steinbeck
 *@created        4. Mai 2005
 *@cdk.created    2002-02-24
 *@cdk.keyword    permutation
 */
public class AtomContainerBondPermutor extends AtomContainerPermutor
{
	
	public AtomContainerBondPermutor(AtomContainer ac)
	{
		setAtomContainer(ac);
		N = atomContainer.getBondCount();
		initBookkeeping();
		initObjectArray();
	}
	
	public void initObjectArray()
	{
		Bond[] bonds = atomContainer.getBonds();
		objects = new Object[atomContainer.getBondCount()];
		for (int f = 0; f < N; f++)
		{
			objects[f] = bonds[f];	
		}
	}
	
	AtomContainer makeResult()
	{
		Bond[] bonds = new Bond[objects.length];
		for (int f = 0; f < objects.length; f++)
		{
			bonds[f] = ((Bond)objects[f]);	
		}
		AtomContainer ac = new AtomContainer(atomContainer);
		ac.setElectronContainers(bonds);
		return (AtomContainer)ac.clone();
	}
	
}

