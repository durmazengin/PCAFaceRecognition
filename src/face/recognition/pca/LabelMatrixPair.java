package face.recognition.pca;

/*-----------------------------------------------------------------------------
File         : LabelMatrixPair
 Author      : Engin DURMAZ
 Date        : 23.04.2019
 Description : Class to store label and related image data together
-----------------------------------------------------------------------------*/

import java.awt.image.DataBuffer;
import java.io.File;
import java.io.IOException;

public class LabelMatrixPair 
{
	private Matrix matrix = null;
	private String label = "";
	private double distance = 0;

	public LabelMatrixPair()
	{

	}

	public LabelMatrixPair(Matrix matrixOfFile, String labelOfFile)
	{
		// convert [M, N] matrix to [M*N, 1] vector		
		matrix = matrixOfFile.convertToVector();
		
		label = labelOfFile;
	}

	public static LabelMatrixPair fromFile(File pathToImage)
	{
		String labelOfFile = pathToImage.getParentFile().getName();
		
		// get pixels to matrix of M x N
		Matrix matrixOfFile = getMatrixFromFile(pathToImage);
				
		LabelMatrixPair pair = new LabelMatrixPair(matrixOfFile, labelOfFile);

		return pair;
	}


	public Matrix getMatrix() 
	{
		return matrix;
	}

	public void setMatrix(Matrix matrix) 
	{
		this.matrix = matrix;
	}

	public String getLabel() 
	{
		return label;
	}

	public void setLabel(String label) 
	{
		this.label = label;
	}

	public double getDistance() 
	{
		return distance;
	}

	public void setDistance(double newDistance) 
	{
		this.distance = newDistance;
	}
	
	static Matrix getMatrixFromFile(File srcfile) 
	{
		java.awt.image.BufferedImage image = null;
		Matrix matrixOfFile = null;
		
		try {
			image = javax.imageio.ImageIO.read(srcfile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if(null != image)
		{
			int width = image.getWidth();
			int height = image.getHeight();
					      
			matrixOfFile = new Matrix(width, height);

			DataBuffer imgDataBuffer = image.getRaster().getDataBuffer();
			for (int row = 0; row < height; row++) 
			{
				for (int col = 0; col < width; col++) 
				{
					//matrixOfFile.set(row, col, image.getRGB(col, row));
					int pixel = imgDataBuffer.getElem(row * width + col);
					matrixOfFile.set(row, col, pixel);
				}
			}
		}
		
		return matrixOfFile;
	}
	
	/*
	 *  get folder of image collection
	 *  for example .../training_images/02/1.bmp
	 *  1 - remove '1.bmp' and obtain '.../training_images/02'
	 *  then get string after last / and obtain '02'
	 */
	private static String getParentFolder(String fullPath) 
	{
		String parentFolder = fullPath.replace('\\', '/');
		if(parentFolder.indexOf('/') > 0)
		{
			parentFolder = parentFolder.substring(0, parentFolder.lastIndexOf('/'));//remove file name
			parentFolder = parentFolder.substring(parentFolder.lastIndexOf('/') + 1);//get folder
		}
		return parentFolder;
	}
}
