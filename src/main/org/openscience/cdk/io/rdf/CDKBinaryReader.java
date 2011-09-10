/* Copyright (C) 2011  Egon Willighagen <egonw@users.sf.net>
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.io.rdf;

import java.io.IOException;
import java.io.InputStream;
import java.io.Reader;

import org.openscience.cdk.Molecule;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.DefaultChemObjectReader;
import org.openscience.cdk.io.formats.CDKBinaryFormat;
import org.openscience.cdk.io.formats.CDKOWLFormat;
import org.openscience.cdk.io.formats.IResourceFormat;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.ObjectBuffer;

/**
 * Reads content from a binary CDK serialization.
 *
 * @cdk.module  io
 * @cdk.githash
 */
@TestClass("org.openscience.cdk.io.rdf.CDKBinReaderTest")
public class CDKBinaryReader extends DefaultChemObjectReader {

    private InputStream input;

    /**
     * Creates a new CDKOWLReader sending output to the given Writer.
     *
     * @param input {@link Reader} from which is OWL input is taken.
     */
    public CDKBinaryReader(InputStream input) {
        this.input = input;
    }
    
    /**
     * Creates a new CDKOWLReader with an undefined input.
     */
    public CDKBinaryReader() {
        this.input = null;
    }
    
    /**
     * Returns the {@link IResourceFormat} for this reader.
     *
     * @return returns a {@link CDKOWLFormat}.
     */
    @TestMethod("testGetFormat")
    public IResourceFormat getFormat() {
        return CDKBinaryFormat.getInstance();
    }

    /**
     * This method must not be used; XML reading requires the use of an
     * {@link InputStream}. Use {@link #setReader(InputStream)} instead.
     *
     * @param reader reader to which should be written.
     * @deprecated 
     */
    @TestMethod("testSetReader_Reader")
    public void setReader(Reader reader) throws CDKException {
        throw new CDKException("Unsupported");
    }

    /** {@inheritDoc} */
    @TestMethod("testSetReader_InputStream")
    public void setReader(InputStream input) throws CDKException {
        this.input = input;
    }

    /** {@inheritDoc} */
	@TestMethod("testAccepts")
    public boolean accepts(Class<? extends IChemObject> classObject) {
		Class[] interfaces = classObject.getInterfaces();
		for (int i=0; i<interfaces.length; i++) {
            if (IAtomContainer.class.equals(interfaces[i])) return true;
			if (IMolecule.class.equals(interfaces[i])) return true;
		}
		Class superClass = classObject.getSuperclass();
	    if (superClass != null) return this.accepts(superClass);
	    return false;
	}

    /** {@inheritDoc} */
    public <T extends IChemObject> T read(T object) throws CDKException {
        if (!(object instanceof IMolecule))
            throw new CDKException(
                "Only supported is reading of IMolecule objects."
            );
        Kryo kryo = new Kryo();
        kryo.setRegistrationOptional(true);
        ObjectBuffer buffer = new ObjectBuffer(kryo);
        Object readObject = buffer.readClassAndObject(input);
        return (T)readObject;
    }

    /** {@inheritDoc} */
    @TestMethod("testClose")
    public void close() throws IOException {
        input.close();
    }

}

