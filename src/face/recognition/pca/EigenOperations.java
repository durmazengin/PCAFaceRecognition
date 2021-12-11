package face.recognition.pca;
/*-----------------------------------------------------------------------------
File         : EigenOperations
 Author      : Engin DURMAZ
 Date        : 23.04.2019
 Description : Class to calculate linear operations like eigen calculations
-----------------------------------------------------------------------------*/

public class EigenOperations 
{
	int numOfColumns = 0; // number of columns in matrix
	double [] dValues = null; // array of eigenvalues
	double [] eValues = null; // array of eigenvalues
	double [][] matrixOfEigenVectors = null; // eigenvectors
	
	double [][] matrixOfHessenbergForm = null; // nonsymmetric Hessenberg form
	boolean isSymetric = false;
	
	double [] orthogonals = null;

	/***************************************************************/
	public EigenOperations (Matrix srcMatrix) 
	{
		double[][] A = new double[srcMatrix.getRowCount()][srcMatrix.getColumnCount()];
		for(int i = 0; i < srcMatrix.getRowCount(); i++)
		{
			for(int j = 0; j < srcMatrix.getColumnCount(); j++)
			{
				A[i][j] = srcMatrix.get(i, j);
			}	
		}
		
		numOfColumns = srcMatrix.getColumnCount();
		matrixOfEigenVectors = new double[numOfColumns][numOfColumns];
		dValues = new double[numOfColumns];
		eValues = new double[numOfColumns];

		isSymetric = true;
		for (int j = 0; j < numOfColumns; j++) 
		{
			for (int i = 0; i < numOfColumns; i++) 
			{
				isSymetric = (A[i][j] == A[j][i]);
				if(!isSymetric)
				{
					break;
				}
			}
			if(!isSymetric)
			{
				break;
			}
		}

		if (isSymetric) 
		{
			for (int i = 0; i < numOfColumns; i++) 
			{
				for (int j = 0; j < numOfColumns; j++) 
				{
					matrixOfEigenVectors[i][j] = A[i][j];
				}
			}

			// Tridiagonalize.
			tridiagonalize();

			// Diagonalize.
			diagonalize();

		} 
		else 
		{
			matrixOfHessenbergForm = new double[numOfColumns][numOfColumns];
			orthogonals = new double[numOfColumns];

			for (int j = 0; j < numOfColumns; j++) 
			{
				for (int i = 0; i < numOfColumns; i++) 
				{
					matrixOfHessenbergForm[i][j] = A[i][j];
				}
			}

			// Reduce to Hessenberg form.
			findOrthogonals();

			// Reduce Hessenberg to real Schur form.
			reduceToSchurForm();
		}
	}

	public double[] getd() 
	{
		return dValues;
	}

	public Matrix getV() 
	{
		Matrix matrixV = new Matrix(numOfColumns, numOfColumns);
		for(int i = 0; i < numOfColumns; i++)
		{
			for(int j = 0; j < numOfColumns; j++)
			{
				matrixV.set(i, j, matrixOfEigenVectors[i][j]);
			}
		}
		return matrixV;
	}
	
	/********************** OUTSOURCE FUNCTIONS ****************/
	private void tridiagonalize() 
	{
		/*
		 *  Source :  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
		 *  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
		 *  Fortran subroutine in EISPACK.
		 */

		for (int j = 0; j < numOfColumns; j++) 
		{
			dValues[j] = matrixOfEigenVectors[numOfColumns-1][j];
		}

		// Householder reduction to tridiagonal form.

		for (int i = numOfColumns-1; i > 0; i--) {

			// Scale to avoid under/overflow.

			double scale = 0.0;
			double h = 0.0;
			for (int k = 0; k < i; k++) {
				scale = scale + Math.abs(dValues[k]);
			}
			if (scale == 0.0) {
				eValues[i] = dValues[i-1];
				for (int j = 0; j < i; j++) {
					dValues[j] = matrixOfEigenVectors[i-1][j];
					matrixOfEigenVectors[i][j] = 0.0;
					matrixOfEigenVectors[j][i] = 0.0;
				}
			} else {

				// Generate Householder vector.

				for (int k = 0; k < i; k++) {
					dValues[k] /= scale;
					h += dValues[k] * dValues[k];
				}
				double f = dValues[i-1];
				double g = Math.sqrt(h);
				if (f > 0) {
					g = -g;
				}
				eValues[i] = scale * g;
				h = h - f * g;
				dValues[i-1] = f - g;
				for (int j = 0; j < i; j++) {
					eValues[j] = 0.0;
				}

				// Apply similarity transformation to remaining columns.

				for (int j = 0; j < i; j++) {
					f = dValues[j];
					matrixOfEigenVectors[j][i] = f;
					g = eValues[j] + matrixOfEigenVectors[j][j] * f;
					for (int k = j+1; k <= i-1; k++) {
						g += matrixOfEigenVectors[k][j] * dValues[k];
						eValues[k] += matrixOfEigenVectors[k][j] * f;
					}
					eValues[j] = g;
				}
				f = 0.0;
				for (int j = 0; j < i; j++) {
					eValues[j] /= h;
					f += eValues[j] * dValues[j];
				}
				double hh = f / (h + h);
				for (int j = 0; j < i; j++) {
					eValues[j] -= hh * dValues[j];
				}
				for (int j = 0; j < i; j++) {
					f = dValues[j];
					g = eValues[j];
					for (int k = j; k <= i-1; k++) {
						matrixOfEigenVectors[k][j] -= (f * eValues[k] + g * dValues[k]);
					}
					dValues[j] = matrixOfEigenVectors[i-1][j];
					matrixOfEigenVectors[i][j] = 0.0;
				}
			}
			dValues[i] = h;
		}

		// Accumulate transformations.

		for (int i = 0; i < numOfColumns-1; i++) {
			matrixOfEigenVectors[numOfColumns-1][i] = matrixOfEigenVectors[i][i];
			matrixOfEigenVectors[i][i] = 1.0;
			double h = dValues[i+1];
			if (h != 0.0) {
				for (int k = 0; k <= i; k++) {
					dValues[k] = matrixOfEigenVectors[k][i+1] / h;
				}
				for (int j = 0; j <= i; j++) {
					double g = 0.0;
					for (int k = 0; k <= i; k++) {
						g += matrixOfEigenVectors[k][i+1] * matrixOfEigenVectors[k][j];
					}
					for (int k = 0; k <= i; k++) {
						matrixOfEigenVectors[k][j] -= g * dValues[k];
					}
				}
			}
			for (int k = 0; k <= i; k++) {
				matrixOfEigenVectors[k][i+1] = 0.0;
			}
		}
		for (int j = 0; j < numOfColumns; j++) {
			dValues[j] = matrixOfEigenVectors[numOfColumns-1][j];
			matrixOfEigenVectors[numOfColumns-1][j] = 0.0;
		}
		matrixOfEigenVectors[numOfColumns-1][numOfColumns-1] = 1.0;
		eValues[0] = 0.0;
	} 

	// diagonalize.

	private void diagonalize() 
	{

		/*
		 *  Source :  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
		 *  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
		 *  Fortran subroutine in EISPACK.
		 */

		for (int i = 1; i < numOfColumns; i++) 
		{
			eValues[i-1] = eValues[i];
		}
		eValues[numOfColumns-1] = 0.0;

		double f = 0.0;
		double tst1 = 0.0;
		double eps = Math.pow(2.0,-52.0);
		for (int l = 0; l < numOfColumns; l++) {

			// Find small subdiagonal element

			tst1 = Math.max(tst1,Math.abs(dValues[l]) + Math.abs(eValues[l]));
			int m = l;
			while (m < numOfColumns) {
				if (Math.abs(eValues[m]) <= eps*tst1) {
					break;
				}
				m++;
			}

			// If m == l, d[l] is an eigenvalue,
					// otherwise, iterate.

			if (m > l) {
				int iter = 0;
				do {
					iter = iter + 1;  // (Could check iteration count here.)

					// Compute implicit shift

					double g = dValues[l];
					double p = (dValues[l+1] - g) / (2.0 * eValues[l]);
					double r = Math.hypot(p, 1.0); // if error check this
					if (p < 0) {
						r = -r;
					}
					dValues[l] = eValues[l] / (p + r);
					dValues[l+1] = eValues[l] * (p + r);
					double dl1 = dValues[l+1];
					double h = g - dValues[l];
					for (int i = l+2; i < numOfColumns; i++) {
						dValues[i] -= h;
					}
					f = f + h;

					// Implicit QL transformation.

					p = dValues[m];
					double c = 1.0;
					double c2 = c;
					double c3 = c;
					double el1 = eValues[l+1];
					double s = 0.0;
					double s2 = 0.0;
					for (int i = m-1; i >= l; i--) {
						c3 = c2;
						c2 = c;
						s2 = s;
						g = c * eValues[i];
						h = c * p;
						r = Math.hypot(p,eValues[i]);// if error check this
						eValues[i+1] = s * r;
						s = eValues[i] / r;
						c = p / r;
						p = c * dValues[i] - s * g;
						dValues[i+1] = h + s * (c * g + s * dValues[i]);

						// Accumulate transformation.

						for (int k = 0; k < numOfColumns; k++) {
							h = matrixOfEigenVectors[k][i+1];
							matrixOfEigenVectors[k][i+1] = s * matrixOfEigenVectors[k][i] + c * h;
							matrixOfEigenVectors[k][i] = c * matrixOfEigenVectors[k][i] - s * h;
						}
					}
					p = -s * s2 * c3 * el1 * eValues[l] / dl1;
					eValues[l] = s * p;
					dValues[l] = c * p;

					// Check for convergence.

				} while (Math.abs(eValues[l]) > eps*tst1);
			}
			dValues[l] = dValues[l] + f;
			eValues[l] = 0.0;
		}

		// Sort eigenvalues and corresponding vectors.

		for (int i = 0; i < numOfColumns-1; i++) {
			int k = i;
			double p = dValues[i];
			for (int j = i+1; j < numOfColumns; j++) {
				if (dValues[j] < p) {
					k = j;
					p = dValues[j];
				}
			}
			if (k != i) {
				dValues[k] = dValues[i];
				dValues[i] = p;
				for (int j = 0; j < numOfColumns; j++) {
					p = matrixOfEigenVectors[j][i];
					matrixOfEigenVectors[j][i] = matrixOfEigenVectors[j][k];
					matrixOfEigenVectors[j][k] = p;
				}
			}
		}
	}

	// Nonsymmetric reduction to Hessenberg form.

	private void findOrthogonals () 
	{
		/*
		 *  Source :  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
		 *  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
		 *  Fortran subroutine in EISPACK.
		 */

		int low = 0;
		int high = numOfColumns-1;

		for (int m = low+1; m <= high-1; m++) {

			// Scale column.

			double scale = 0.0;
			for (int i = m; i <= high; i++) {
				scale = scale + Math.abs(matrixOfHessenbergForm[i][m-1]);
			}
			if (scale != 0.0) {

				// Compute Householder transformation.

				double h = 0.0;
				for (int i = high; i >= m; i--) {
					orthogonals[i] = matrixOfHessenbergForm[i][m-1]/scale;
					h += orthogonals[i] * orthogonals[i];
				}
				double g = Math.sqrt(h);
				if (orthogonals[m] > 0) {
					g = -g;
				}
				h = h - orthogonals[m] * g;
				orthogonals[m] = orthogonals[m] - g;

				// Apply Householder similarity transformation
				// H = (I-u*u'/h)*H*(I-u*u')/h)

				for (int j = m; j < numOfColumns; j++) {
					double f = 0.0;
					for (int i = high; i >= m; i--) {
						f += orthogonals[i]*matrixOfHessenbergForm[i][j];
					}
					f = f/h;
					for (int i = m; i <= high; i++) {
						matrixOfHessenbergForm[i][j] -= f*orthogonals[i];
					}
				}

				for (int i = 0; i <= high; i++) {
					double f = 0.0;
					for (int j = high; j >= m; j--) {
						f += orthogonals[j]*matrixOfHessenbergForm[i][j];
					}
					f = f/h;
					for (int j = m; j <= high; j++) {
						matrixOfHessenbergForm[i][j] -= f*orthogonals[j];
					}
				}
				orthogonals[m] = scale*orthogonals[m];
				matrixOfHessenbergForm[m][m-1] = scale*g;
			}
		}

		// Accumulate transformations (Algol's ortran).

		for (int i = 0; i < numOfColumns; i++) {
			for (int j = 0; j < numOfColumns; j++) {
				matrixOfEigenVectors[i][j] = (i == j ? 1.0 : 0.0);
			}
		}

		for (int m = high-1; m >= low+1; m--) {
			if (matrixOfHessenbergForm[m][m-1] != 0.0) {
				for (int i = m+1; i <= high; i++) {
					orthogonals[i] = matrixOfHessenbergForm[i][m-1];
				}
				for (int j = m; j <= high; j++) {
					double g = 0.0;
					for (int i = m; i <= high; i++) {
						g += orthogonals[i] * matrixOfEigenVectors[i][j];
					}
					// Double division avoids possible underflow
					g = (g / orthogonals[m]) / matrixOfHessenbergForm[m][m-1];
					for (int i = m; i <= high; i++) {
						matrixOfEigenVectors[i][j] += g * orthogonals[i];
					}
				}
			}
		}
	}


	// Complex scalar division.

	private transient double cdivr, cdivi;
	private void cdiv(double xr, double xi, double yr, double yi) {
		double r,d;
		if (Math.abs(yr) > Math.abs(yi)) {
			r = yi/yr;
			d = yr + r*yi;
			cdivr = (xr + r*xi)/d;
			cdivi = (xi - r*xr)/d;
		} else {
			r = yr/yi;
			d = yi + r*yr;
			cdivr = (r*xr + xi)/d;
			cdivi = (r*xi - xr)/d;
		}
	}


	// Nonsymmetric reduction from Hessenberg to real Schur form.

	private void reduceToSchurForm () 
	{

		/*
		 *  Source :  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
		 *  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
		 *  Fortran subroutine in EISPACK.
		 */

		// Initialize

		int nn = this.numOfColumns;
		int n = nn-1;
		int low = 0;
		int high = nn-1;
		double eps = Math.pow(2.0,-52.0);
		double exshift = 0.0;
		double p=0,q=0,r=0,s=0,z=0,t,w,x,y;

		// Store roots isolated by balanc and compute matrix norm

		double norm = 0.0;
		for (int i = 0; i < nn; i++) {
			if (i < low | i > high) {
				dValues[i] = matrixOfHessenbergForm[i][i];
				eValues[i] = 0.0;
			}
			for (int j = Math.max(i-1,0); j < nn; j++) {
				norm = norm + Math.abs(matrixOfHessenbergForm[i][j]);
			}
		}

		// Outer loop over eigenvalue index

		int iter = 0;
		while (n >= low) {

			// Look for single small sub-diagonal element

			int l = n;
			while (l > low) {
				s = Math.abs(matrixOfHessenbergForm[l-1][l-1]) + Math.abs(matrixOfHessenbergForm[l][l]);
				if (s == 0.0) {
					s = norm;
				}
				if (Math.abs(matrixOfHessenbergForm[l][l-1]) < eps * s) {
					break;
				}
				l--;
			}

			// Check for convergence
			// One root found

			if (l == n) {
				matrixOfHessenbergForm[n][n] = matrixOfHessenbergForm[n][n] + exshift;
				dValues[n] = matrixOfHessenbergForm[n][n];
				eValues[n] = 0.0;
				n--;
				iter = 0;

				// Two roots found

			} else if (l == n-1) {
				w = matrixOfHessenbergForm[n][n-1] * matrixOfHessenbergForm[n-1][n];
				p = (matrixOfHessenbergForm[n-1][n-1] - matrixOfHessenbergForm[n][n]) / 2.0;
				q = p * p + w;
				z = Math.sqrt(Math.abs(q));
				matrixOfHessenbergForm[n][n] = matrixOfHessenbergForm[n][n] + exshift;
				matrixOfHessenbergForm[n-1][n-1] = matrixOfHessenbergForm[n-1][n-1] + exshift;
				x = matrixOfHessenbergForm[n][n];

				// Real pair

				if (q >= 0) {
					if (p >= 0) {
						z = p + z;
					} else {
						z = p - z;
					}
					dValues[n-1] = x + z;
					dValues[n] = dValues[n-1];
					if (z != 0.0) {
						dValues[n] = x - w / z;
					}
					eValues[n-1] = 0.0;
					eValues[n] = 0.0;
					x = matrixOfHessenbergForm[n][n-1];
					s = Math.abs(x) + Math.abs(z);
					p = x / s;
					q = z / s;
					r = Math.sqrt(p * p+q * q);
					p = p / r;
					q = q / r;

					// Row modification

					for (int j = n-1; j < nn; j++) {
						z = matrixOfHessenbergForm[n-1][j];
						matrixOfHessenbergForm[n-1][j] = q * z + p * matrixOfHessenbergForm[n][j];
						matrixOfHessenbergForm[n][j] = q * matrixOfHessenbergForm[n][j] - p * z;
					}

					// Column modification

					for (int i = 0; i <= n; i++) {
						z = matrixOfHessenbergForm[i][n-1];
						matrixOfHessenbergForm[i][n-1] = q * z + p * matrixOfHessenbergForm[i][n];
						matrixOfHessenbergForm[i][n] = q * matrixOfHessenbergForm[i][n] - p * z;
					}

					// Accumulate transformations

					for (int i = low; i <= high; i++) {
						z = matrixOfEigenVectors[i][n-1];
						matrixOfEigenVectors[i][n-1] = q * z + p * matrixOfEigenVectors[i][n];
						matrixOfEigenVectors[i][n] = q * matrixOfEigenVectors[i][n] - p * z;
					}

					// Complex pair

				} else {
					dValues[n-1] = x + p;
					dValues[n] = x + p;
					eValues[n-1] = z;
					eValues[n] = -z;
				}
				n = n - 2;
				iter = 0;

				// No convergence yet

			} else {

				// Form shift

				x = matrixOfHessenbergForm[n][n];
				y = 0.0;
				w = 0.0;
				if (l < n) {
					y = matrixOfHessenbergForm[n-1][n-1];
					w = matrixOfHessenbergForm[n][n-1] * matrixOfHessenbergForm[n-1][n];
				}

				// Wilkinson's original ad hoc shift

				if (iter == 10) {
					exshift += x;
					for (int i = low; i <= n; i++) {
						matrixOfHessenbergForm[i][i] -= x;
					}
					s = Math.abs(matrixOfHessenbergForm[n][n-1]) + Math.abs(matrixOfHessenbergForm[n-1][n-2]);
					x = y = 0.75 * s;
					w = -0.4375 * s * s;
				}

				// MATLAB's new ad hoc shift

				if (iter == 30) {
					s = (y - x) / 2.0;
					s = s * s + w;
					if (s > 0) {
						s = Math.sqrt(s);
						if (y < x) {
							s = -s;
						}
						s = x - w / ((y - x) / 2.0 + s);
						for (int i = low; i <= n; i++) {
							matrixOfHessenbergForm[i][i] -= s;
						}
						exshift += s;
						x = y = w = 0.964;
					}
				}

				iter = iter + 1;   // (Could check iteration count here.)

				// Look for two consecutive small sub-diagonal elements

				int m = n-2;
				while (m >= l) {
					z = matrixOfHessenbergForm[m][m];
					r = x - z;
					s = y - z;
					p = (r * s - w) / matrixOfHessenbergForm[m+1][m] + matrixOfHessenbergForm[m][m+1];
					q = matrixOfHessenbergForm[m+1][m+1] - z - r - s;
					r = matrixOfHessenbergForm[m+2][m+1];
					s = Math.abs(p) + Math.abs(q) + Math.abs(r);
					p = p / s;
					q = q / s;
					r = r / s;
					if (m == l) {
						break;
					}
					if (Math.abs(matrixOfHessenbergForm[m][m-1]) * (Math.abs(q) + Math.abs(r)) <
							eps * (Math.abs(p) * (Math.abs(matrixOfHessenbergForm[m-1][m-1]) + Math.abs(z) +
									Math.abs(matrixOfHessenbergForm[m+1][m+1])))) {
						break;
					}
					m--;
				}

				for (int i = m+2; i <= n; i++) {
					matrixOfHessenbergForm[i][i-2] = 0.0;
					if (i > m+2) {
						matrixOfHessenbergForm[i][i-3] = 0.0;
					}
				}

				// Double QR step involving rows l:n and columns m:n


				for (int k = m; k <= n-1; k++) {
					boolean notlast = (k != n-1);
					if (k != m) {
						p = matrixOfHessenbergForm[k][k-1];
						q = matrixOfHessenbergForm[k+1][k-1];
						r = (notlast ? matrixOfHessenbergForm[k+2][k-1] : 0.0);
						x = Math.abs(p) + Math.abs(q) + Math.abs(r);
						if (x == 0.0) {
							continue;
						}
						p = p / x;
						q = q / x;
						r = r / x;
					}

					s = Math.sqrt(p * p + q * q + r * r);
					if (p < 0) {
						s = -s;
					}
					if (s != 0) {
						if (k != m) {
							matrixOfHessenbergForm[k][k-1] = -s * x;
						} else if (l != m) {
							matrixOfHessenbergForm[k][k-1] = -matrixOfHessenbergForm[k][k-1];
						}
						p = p + s;
						x = p / s;
						y = q / s;
						z = r / s;
						q = q / p;
						r = r / p;

						// Row modification

						for (int j = k; j < nn; j++) {
							p = matrixOfHessenbergForm[k][j] + q * matrixOfHessenbergForm[k+1][j];
							if (notlast) {
								p = p + r * matrixOfHessenbergForm[k+2][j];
								matrixOfHessenbergForm[k+2][j] = matrixOfHessenbergForm[k+2][j] - p * z;
							}
							matrixOfHessenbergForm[k][j] = matrixOfHessenbergForm[k][j] - p * x;
							matrixOfHessenbergForm[k+1][j] = matrixOfHessenbergForm[k+1][j] - p * y;
						}

						// Column modification

						for (int i = 0; i <= Math.min(n,k+3); i++) {
							p = x * matrixOfHessenbergForm[i][k] + y * matrixOfHessenbergForm[i][k+1];
							if (notlast) {
								p = p + z * matrixOfHessenbergForm[i][k+2];
								matrixOfHessenbergForm[i][k+2] = matrixOfHessenbergForm[i][k+2] - p * r;
							}
							matrixOfHessenbergForm[i][k] = matrixOfHessenbergForm[i][k] - p;
							matrixOfHessenbergForm[i][k+1] = matrixOfHessenbergForm[i][k+1] - p * q;
						}

						// Accumulate transformations

						for (int i = low; i <= high; i++) {
							p = x * matrixOfEigenVectors[i][k] + y * matrixOfEigenVectors[i][k+1];
							if (notlast) {
								p = p + z * matrixOfEigenVectors[i][k+2];
								matrixOfEigenVectors[i][k+2] = matrixOfEigenVectors[i][k+2] - p * r;
							}
							matrixOfEigenVectors[i][k] = matrixOfEigenVectors[i][k] - p;
							matrixOfEigenVectors[i][k+1] = matrixOfEigenVectors[i][k+1] - p * q;
						}
					}  // (s != 0)
				}  // k loop
			}  // check convergence
		}  // while (n >= low)

		// Backsubstitute to find vectors of upper triangular form

		if (norm == 0.0) {
			return;
		}

		for (n = nn-1; n >= 0; n--) {
			p = dValues[n];
			q = eValues[n];

			// Real vector

			if (q == 0) {
				int l = n;
				matrixOfHessenbergForm[n][n] = 1.0;
				for (int i = n-1; i >= 0; i--) {
					w = matrixOfHessenbergForm[i][i] - p;
					r = 0.0;
					for (int j = l; j <= n; j++) {
						r = r + matrixOfHessenbergForm[i][j] * matrixOfHessenbergForm[j][n];
					}
					if (eValues[i] < 0.0) {
						z = w;
						s = r;
					} else {
						l = i;
						if (eValues[i] == 0.0) {
							if (w != 0.0) {
								matrixOfHessenbergForm[i][n] = -r / w;
							} else {
								matrixOfHessenbergForm[i][n] = -r / (eps * norm);
							}

							// Solve real equations

						} else {
							x = matrixOfHessenbergForm[i][i+1];
							y = matrixOfHessenbergForm[i+1][i];
							q = (dValues[i] - p) * (dValues[i] - p) + eValues[i] * eValues[i];
							t = (x * s - z * r) / q;
							matrixOfHessenbergForm[i][n] = t;
							if (Math.abs(x) > Math.abs(z)) {
								matrixOfHessenbergForm[i+1][n] = (-r - w * t) / x;
							} else {
								matrixOfHessenbergForm[i+1][n] = (-s - y * t) / z;
							}
						}

						// Overflow control

						t = Math.abs(matrixOfHessenbergForm[i][n]);
						if ((eps * t) * t > 1) {
							for (int j = i; j <= n; j++) {
								matrixOfHessenbergForm[j][n] = matrixOfHessenbergForm[j][n] / t;
							}
						}
					}
				}

				// Complex vector

			} else if (q < 0) {
				int l = n-1;

				// Last vector component imaginary so matrix is triangular

				if (Math.abs(matrixOfHessenbergForm[n][n-1]) > Math.abs(matrixOfHessenbergForm[n-1][n])) {
					matrixOfHessenbergForm[n-1][n-1] = q / matrixOfHessenbergForm[n][n-1];
					matrixOfHessenbergForm[n-1][n] = -(matrixOfHessenbergForm[n][n] - p) / matrixOfHessenbergForm[n][n-1];
				} else {
					cdiv(0.0,-matrixOfHessenbergForm[n-1][n],matrixOfHessenbergForm[n-1][n-1]-p,q);
					matrixOfHessenbergForm[n-1][n-1] = cdivr;
					matrixOfHessenbergForm[n-1][n] = cdivi;
				}
				matrixOfHessenbergForm[n][n-1] = 0.0;
				matrixOfHessenbergForm[n][n] = 1.0;
				for (int i = n-2; i >= 0; i--) {
					double ra,sa,vr,vi;
					ra = 0.0;
					sa = 0.0;
					for (int j = l; j <= n; j++) {
						ra = ra + matrixOfHessenbergForm[i][j] * matrixOfHessenbergForm[j][n-1];
						sa = sa + matrixOfHessenbergForm[i][j] * matrixOfHessenbergForm[j][n];
					}
					w = matrixOfHessenbergForm[i][i] - p;

					if (eValues[i] < 0.0) {
						z = w;
						r = ra;
						s = sa;
					} else {
						l = i;
						if (eValues[i] == 0) {
							cdiv(-ra,-sa,w,q);
							matrixOfHessenbergForm[i][n-1] = cdivr;
							matrixOfHessenbergForm[i][n] = cdivi;
						} else {

							// Solve complex equations

							x = matrixOfHessenbergForm[i][i+1];
							y = matrixOfHessenbergForm[i+1][i];
							vr = (dValues[i] - p) * (dValues[i] - p) + eValues[i] * eValues[i] - q * q;
							vi = (dValues[i] - p) * 2.0 * q;
							if (vr == 0.0 & vi == 0.0) {
								vr = eps * norm * (Math.abs(w) + Math.abs(q) +
										Math.abs(x) + Math.abs(y) + Math.abs(z));
							}
							cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi);
							matrixOfHessenbergForm[i][n-1] = cdivr;
							matrixOfHessenbergForm[i][n] = cdivi;
							if (Math.abs(x) > (Math.abs(z) + Math.abs(q))) {
								matrixOfHessenbergForm[i+1][n-1] = (-ra - w * matrixOfHessenbergForm[i][n-1] + q * matrixOfHessenbergForm[i][n]) / x;
								matrixOfHessenbergForm[i+1][n] = (-sa - w * matrixOfHessenbergForm[i][n] - q * matrixOfHessenbergForm[i][n-1]) / x;
							} else {
								cdiv(-r-y*matrixOfHessenbergForm[i][n-1],-s-y*matrixOfHessenbergForm[i][n],z,q);
								matrixOfHessenbergForm[i+1][n-1] = cdivr;
								matrixOfHessenbergForm[i+1][n] = cdivi;
							}
						}

						// Overflow control

						t = Math.max(Math.abs(matrixOfHessenbergForm[i][n-1]),Math.abs(matrixOfHessenbergForm[i][n]));
						if ((eps * t) * t > 1) {
							for (int j = i; j <= n; j++) {
								matrixOfHessenbergForm[j][n-1] = matrixOfHessenbergForm[j][n-1] / t;
								matrixOfHessenbergForm[j][n] = matrixOfHessenbergForm[j][n] / t;
							}
						}
					}
				}
			}
		}

		// Vectors of isolated roots

		for (int i = 0; i < nn; i++) {
			if (i < low | i > high) {
				for (int j = i; j < nn; j++) {
					matrixOfEigenVectors[i][j] = matrixOfHessenbergForm[i][j];
				}
			}
		}

		// Back transformation to get eigenvectors of original matrix

		for (int j = nn-1; j >= low; j--) {
			for (int i = low; i <= high; i++) {
				z = 0.0;
				for (int k = low; k <= Math.min(j,high); k++) {
					z = z + matrixOfEigenVectors[i][k] * matrixOfHessenbergForm[k][j];
				}
				matrixOfEigenVectors[i][j] = z;
			}
		}
	}

}
