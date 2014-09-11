/* Copyright (C) 2005-2008  Nina Jeliazkova <nina@acad.bg>
 * 
 * Contact: cdk-devel@lists.sourceforge.net
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package org.openscience.cdk.qsar.descriptors.molecular;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;

import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.IsomorphismTester;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.qsar.AtomValenceTool;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.cdk.smiles.smarts.parser.SMARTSParser;
import org.openscience.cdk.templates.MoleculeFactory;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class PKASmartsDescriptorTest {

	private PKASmartsDescriptor pka;
	private final CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
	private SmilesParser smiParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
	private Aromaticity cdk = new Aromaticity(ElectronDonation.cdk(), Cycles.all());

    @Before
    public void setUp() throws Exception {
    	 pka = new PKASmartsDescriptor();
    }

    @After
    public void tearDown() throws Exception {
    	pka = null;
    }

    @Test
    public void test() throws Exception {
        Hashtable<Integer,PKANode> tree = pka.getTree();
        Enumeration<Integer> k = tree.keys();
        ArrayList<String> smarts = new ArrayList<String>();
        int failedNodes = 0;
        int failedSmarts = 0;
        int nullSmarts = 0;
        int allnodes =0;
        while (k.hasMoreElements()) {
            PKANode node = tree.get(k.nextElement());
            allnodes++;
            try {
                if (node.getSmarts() == null) {
                    nullSmarts++;
                } else {
                	new SMARTSQueryTool(node.getSmarts(), SilentChemObjectBuilder.getInstance());
                }
            } catch (IllegalArgumentException x) {
                failedNodes ++;
                if (smarts.indexOf(node.getSmarts())<0) {
                    smarts.add(node.getSmarts());
                    failedSmarts++;
                }
            }
        }
        Assert.assertEquals(1,nullSmarts); //root smarts
        Assert.assertEquals(0, failedSmarts);
        Assert.assertEquals(0, failedNodes);
        Assert.assertEquals(1527,allnodes);
        
    }

    @Test
    public void testAcidPkaNode() throws Exception {
    	 IAtomContainer mol = smiParser.parseSmiles("O[N+](=O)[O-]");
    	 AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
    	 hAdder.addImplicitHydrogens(mol);
    	 AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
    	 cdk.apply(mol);
    	 for (int i=0; i < mol.getAtomCount();i++) {
    		 IAtom atom = mol.getAtom(i);
    		 atom.setValency(AtomValenceTool.getValence(atom));
    	 }

    	 PKANode node = new PKANode();
    	 node.setSmarts("[G14;H][^1,^2]");
    	 Assert.assertTrue(node.find(mol));
    }
    
    @Test
    public void test462() throws Exception {
    	IAtomContainer mol = smiParser.parseSmiles("O[N+](=O)[O-]");
    	AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
    	hAdder.addImplicitHydrogens(mol);
    	AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
    	cdk.apply(mol);
    	for (int i=0; i < mol.getAtomCount();i++) {
    		IAtom atom = mol.getAtom(i);
    		atom.setValency(AtomValenceTool.getValence(atom));
    	}

    	PKANode node = new PKANode();
    	node.setSmarts("[^1,^2][G14v2]");
    	Assert.assertFalse(node.find(mol));
    }

    @Test
    public void testAcid15()  throws Exception {
    	String smarts = "[^1,^2][G14v2]";
    	String smiles = "O[N+](=O)[O-]";

    	IQueryAtomContainer query = SMARTSParser.parse(smarts, SilentChemObjectBuilder.getInstance());
    	
    	IAtomContainer mol = smiParser.parseSmiles(smiles);
    	AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
    	hAdder.addImplicitHydrogens(mol);
    	AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);

    	IsomorphismTester isoTester = new IsomorphismTester();
    	Assert.assertFalse(isoTester.isIsomorphic(mol, query));
    }

    @Test
    public void testAcid() throws Exception {
    	IAtomContainer mol = smiParser.parseSmiles("O[N+](=O)[O-]");
    	AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
    	hAdder.addImplicitHydrogens(mol);
    	AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
    	cdk.apply(mol);

		for (int i=0; i < mol.getAtomCount();i++) {
    		IAtom atom = mol.getAtom(i);
        	atom.setValency(AtomValenceTool.getValence(atom));
    		Assert.assertNotNull(atom.getValency());
		}
		for (int i=0; i < mol.getBondCount();i++) {
			Assert.assertNotNull(mol.getBond(i).getOrder());
		}
		DescriptorValue value = pka.calculate(mol);
		Assert.assertNotNull(value.getValue());
		// the next doesn't work with the CDK yet (but hopefully in the future!)
		// Assert.assertNotNull(((DescriptorResult)value.getValue()).getExplanation());
    }    

    @Test
    public void testOne() throws Exception {
    	IAtomContainer mol = MoleculeFactory.makeBenzene();
    	AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
    	hAdder.addImplicitHydrogens(mol);
    	AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
    	cdk.apply(mol);
    	
    	for (int i=0; i < mol.getAtomCount();i++) {
    		IAtom atom = mol.getAtom(i);
        	atom.setValency(AtomValenceTool.getValence(atom));
    		Assert.assertNotNull(atom.getValency());
    	}
    	for (int i=0; i < mol.getBondCount();i++) {
    		Assert.assertNotNull(mol.getBond(i).getOrder());
    	}
    	DescriptorValue value = pka.calculate(mol);
    	Assert.assertNotNull(value.getValue());
    }

}
