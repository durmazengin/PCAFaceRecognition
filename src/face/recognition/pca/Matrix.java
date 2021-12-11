package face.recognition.pca;

/*-----------------------------------------------------------------------------
File         : Matrix
 Author      : Engin DURMAZ
 Date        : 23.04.2019
 Description : Class to store matrix data and some matrix operations
              - put a value to matrix, get a value from matrix
              - get row count, get column count
              - calculate Euclidean distance to another matrix
              - transpose operation
              - sum to another matrix
              - subtract another matrix from this
              - multiply with another matrix
-----------------------------------------------------------------------------*/

public class Matrix 
{
	private double[][] values;
	private int numOfRows = 0;
	private int numOfColumns = 0;

	/*
	 * Constructor to create MxN matrix
	 * M = rowCount
	 * N = columnCount
	 */
	public Matrix (int rowCount, int columnCount) 
	{
		this.numOfRows = rowCount;
		this.numOfColumns = columnCount;
		values = new double[rowCount][columnCount];
	}

	/*
	 * Get all values of matrix
	 */
	public double[][] getInArray () {
		return values;
	}

	/*
	 * Get M = number of rows
	 */
	public int getRowCount () {
		return numOfRows;
	}

	/*
	 * Get N = number of columns
	 */
	public int getColumnCount () {
		return numOfColumns;
	}

	/*
	 * Get a value from matrix
	 */
	public double get (int i, int j) 
	{
		return values[i][j];
	}

	/*
	 * Set a value in matrix
	 */
	public boolean set (int i, int j, double newValue) 
	{
		if(i < numOfRows && j < numOfColumns)
		{
			values[i][j] = newValue;
			return true;
		}
		return false;
	}

	/*
	 * Function to calculate Euclidean distance to another matrix
	 */
	public double getEuclideanDistance(Matrix target) {
		int size = this.getRowCount();
		double sum = 0;

		for (int i = 0; i < size; i++) 
		{
			sum += Math.pow(this.get(i, 0) - target.get(i, 0), 2);
		}

		return Math.sqrt(sum);
	}

	/*
	 * Return transpose of this matrix 
	 * Not affect current matrix
	 */
	public Matrix transpose () 
	{		
		Matrix transposeMatrix = new Matrix(numOfColumns, numOfRows);
		
		for (int i = 0; i < numOfRows; i++) 
		{
			for (int j = 0; j < numOfColumns; j++) 
			{
				transposeMatrix.set(j, i, this.get(i, j));
			}
		}
		return transposeMatrix;
	}

	/*
	 * Return sum of given matrix and this matrix
	 * Not affect current matrix
	 */
	public Matrix sum(Matrix matrixToAdd) 
	{
		Matrix result = new Matrix(numOfRows, numOfColumns);
		for(int i = 0; i < numOfRows; i++)
		{
			for(int j = 0; j < numOfColumns; j++)
			{
				double value = this.get(i, j) + matrixToAdd.get(i, j);
				result.set(i, j, value);
			}
		}
		return result;
	}

	/*
	 * Return difference of this matrix and given matrix
	 * Not affect current matrix
	 */
	public Matrix subtract(Matrix subtractor) 
	{
		Matrix result = new Matrix(numOfRows, numOfColumns);
		for(int i = 0; i < numOfRows; i++)
		{
			for(int j = 0; j < numOfColumns; j++)
			{
				double value = this.get(i, j) - subtractor.get(i, j);
				result.set(i, j, value);
			}
		}
		return result;
	}
	
	/*
	 * Return multiplication of given matrix and this
	 * Not affect current matrix
	 */
	public Matrix multiply (Matrix multiplier) 
	{
		if (numOfColumns != multiplier.numOfRows) 
		{
			return null; // 
		}
		Matrix result = new Matrix(numOfRows, multiplier.numOfColumns);
		
		/*
		 * 2 5 1 6  x 4 8 2 = 8+15+7+48   16+45+5+24  ..= 78  90 78
		 * 4 2 7 3    3 9 4   16+6+49+24  32+18+35+12 ..  95  97 82
		 *            7 5 6
		 *            8 4 8
		 */
		for (int row = 0; row < numOfRows; row++) 
		{

			for (int clmnMultiplier = 0; clmnMultiplier < multiplier.numOfColumns; clmnMultiplier++) 
			{
				// [2 3 4 5] * [1 4 8 9]^T
				double value = 0;
				for (int clmn = 0; clmn < numOfColumns; clmn++) 
				{
					value += this.get(row, clmn) * multiplier.get(clmn, clmnMultiplier);
				}
				result.set(row, clmnMultiplier, value);
			}
		}
		return result;
	}
	/*
	 * Copy some range from given matrix to this matrix
	 */
	public void copyFrom(int rowStart, int rowEnd, 
			int clmnStart, int clmnEnd, Matrix srcMatirc) 
	{
		for (int i = rowStart; i <= rowEnd; i++) 
		{
            for (int j = clmnStart; j <= clmnEnd; j++) 
            {
            	this.set(i, j, srcMatirc.get(i - rowStart, j - clmnStart));
            }
         }
		
	}
	/*
	 * convert M x N matrix to (M*N x 1) 
	 */
	public Matrix convertToVector() 
	{
        Matrix vectorMatrix = new Matrix(numOfRows * numOfColumns, 1);
        
        for (int column = 0; column < numOfColumns; column++) 
        {
            for (int row = 0; row < numOfRows; row++) 
            {
                vectorMatrix.set(column * numOfRows + row, 0, this.get(row, column));
            }
        }
        return vectorMatrix;
    }
	/*
	 * convert M * N x 1 matrix to (M x N) 
	 */

	public Matrix convertFromVector(int targetColumnCount) 
	{
        int targetRowCount = numOfRows / targetColumnCount;

        Matrix newMatrix = new Matrix(targetRowCount, targetColumnCount);
        
        for (int column = 0; column < targetColumnCount; column++) 
        {
            for (int row = 0; row < targetRowCount; row++) 
            {
            	newMatrix.set(row, column, this.get(column * numOfRows + row, 0));
            }
        }
        return newMatrix;
    }

	/*
	 * this is a test function
	 * i added it to check if multiplication works
	 */
	public static boolean testMultiplication()
	{
		Matrix a = new Matrix(2,  4);
		Matrix b = new Matrix(4,  3);
		
		a.set(0, 0, 2);
		a.set(0, 1, 5);
		a.set(0, 2, 1);
		a.set(0, 3, 6);

		a.set(1, 0, 4);
		a.set(1, 1, 2);
		a.set(1, 2, 7);
		a.set(1, 3, 3);

		b.set(0, 0, 4);
		b.set(0, 1, 8);
		b.set(0, 2, 2);

		b.set(1, 0, 3);
		b.set(1, 1, 9);
		b.set(1, 2, 4);

		b.set(2, 0, 7);
		b.set(2, 1, 5);
		b.set(2, 2, 6);

		b.set(3, 0, 8);
		b.set(3, 1, 4);
		b.set(3, 2, 8);
		
		Matrix result = a.multiply(b);
		if((result.get(0, 0) != 78) ||
		   (result.get(0, 1) != 90) ||
		   (result.get(0, 2) != 78) ||
		   (result.get(1, 0) != 95) ||
		   (result.get(1, 1) != 97) ||
		   (result.get(1, 2) != 82))
		{
			return false;
		}
		return true;
	}

}
