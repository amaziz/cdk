/* Copyright (C) 2011,2013  Egon Willighagen <egonw@users.sf.net>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.renderer.visitor;

import java.awt.Color;
import java.awt.geom.AffineTransform;
import java.io.ByteArrayInputStream;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.renderer.AtomContainerRenderer;
import org.openscience.cdk.renderer.RendererModel;
import org.openscience.cdk.renderer.elements.TextElement;
import org.openscience.cdk.renderer.font.AWTFontManager;
import org.openscience.cdk.renderer.generators.BasicAtomGenerator;
import org.openscience.cdk.renderer.generators.BasicBondGenerator;
import org.openscience.cdk.renderer.generators.BasicSceneGenerator;
import org.openscience.cdk.renderer.generators.IGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.w3c.dom.Document;

/**
 * @cdk.module  test-rendersvg
 * @cdk.githash
 *
 */
public class SVGGeneratorTest {

	private IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
	
	private StructureDiagramGenerator sdg = new StructureDiagramGenerator();
	
	public IMolecule layout(IMolecule molecule) {
		sdg.setMolecule(molecule);
		try {
			sdg.generateCoordinates();
		} catch (Exception e) {
			System.err.println(e);
		}
		return sdg.getMolecule();
	}

	@Test
	public void testEmptyModel() throws Exception {
		IMolecule dummy = builder.newInstance(IMolecule.class);
        SVGGenerator svgGenerator = new SVGGenerator();
		List<IGenerator<IAtomContainer>> generators =
				new ArrayList<IGenerator<IAtomContainer>>();
		generators.add(new BasicSceneGenerator());
		generators.add(new BasicBondGenerator());
		BasicAtomGenerator atomGenerator = new BasicAtomGenerator();
		generators.add(atomGenerator);
		
		AtomContainerRenderer renderer = new AtomContainerRenderer(generators, new AWTFontManager());
        renderer.paint(dummy, svgGenerator);
        String svg = svgGenerator.getResult();
		Assert.assertNotSame(0, svg.length());
		
		 // parse an XML document into a DOM tree
	    DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
	    dbf.setValidating(true);
	    DocumentBuilder parser = dbf.newDocumentBuilder();
	    Document document = parser.parse (new ByteArrayInputStream (svg.getBytes()));
	    Assert.assertNotNull(document);
	}

	@Test
	public void testConstructor() {
		IDrawVisitor visitor = new SVGGenerator();
		Assert.assertNotNull(visitor);
	}

	@Test
	public void testSetFontManager() {
		IDrawVisitor visitor = new SVGGenerator();
		visitor.setFontManager(new AWTFontManager());
		// at least we now know it did not crash...
		Assert.assertNotNull(visitor);
	}

	@Test
	public void testSetRendererModel() {
		IDrawVisitor visitor = new SVGGenerator();
		visitor.setRendererModel(new RendererModel());
		// at least we now know it did not crash...
		Assert.assertNotNull(visitor);
	}

    @Test
	public void testVisit() throws Exception {
    	SVGGenerator svgGenerator = new SVGGenerator();
		svgGenerator.setFontManager(new AWTFontManager());
		svgGenerator.setTransform(new AffineTransform());
		svgGenerator.visit(new TextElement(2, 3, "Foo", Color.BLACK));
		// at least we now know it did not crash...
		Assert.assertNotNull(svgGenerator);
		String svg = svgGenerator.getResult();
		Assert.assertNotSame(0, svg.length());
		System.out.println(svg);
		
		 // parse an XML document into a DOM tree
	    DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
	    DocumentBuilder parser = dbf.newDocumentBuilder();
	    Document document = parser.parse (new ByteArrayInputStream (svg.getBytes()));
	    Assert.assertNotNull(document);
	}

}
