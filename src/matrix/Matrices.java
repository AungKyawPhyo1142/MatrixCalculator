package Matrix;

import javax.swing.*;
import java.awt.*;

public class Matrices {

public static String s;

public static float[][] InverseMatrix (float[][] matrix) {

		float[][] auxiliaryMatrix, invertedMatrix;
		int[] index;

		auxiliaryMatrix = new float[matrix.length][matrix.length];
		invertedMatrix = new float[matrix.length][matrix.length];
		index = new int[matrix.length];

		for (int i = 0; i < matrix.length; ++i) {
			auxiliaryMatrix[i][i] = 1;
		}

		UpperTriangle (matrix, index);

		for (int i = 0; i < (matrix.length - 1); ++i) {
			for (int j = (i + 1); j < matrix.length; ++j) {
				for (int k = 0; k < matrix.length; ++k) {
					auxiliaryMatrix[index[j]][k] -= matrix[index[j]][i] * auxiliaryMatrix[index[i]][k];
				}
			}
		}

		for (int i = 0; i < matrix.length; ++i) {
			invertedMatrix[matrix.length - 1][i] = (auxiliaryMatrix[index[matrix.length - 1]][i] / matrix[index[matrix.length - 1]][matrix.length - 1]);

			for (int j = (matrix.length - 2); j >= 0; --j) {
				invertedMatrix[j][i] = auxiliaryMatrix[index[j]][i];

				for (int k = (j + 1); k < matrix.length; ++k) {
					invertedMatrix[j][i] -= (matrix[index[j]][k] * invertedMatrix[k][i]);
				}

				invertedMatrix[j][i] /= matrix[index[j]][j];
			}
		}

		return (invertedMatrix);
	}

public static void UpperTriangle (float[][] matrix, int[] index) {

		float[] c;
		float c0, c1, pi0, pi1, pj;
		int itmp, k;

		c = new float[matrix.length];

		for (int i = 0; i < matrix.length; ++i) {
			index[i] = i;
		}

		for (int i = 0; i < matrix.length; ++i) {
			c1 = 0;

			for (int j = 0; j < matrix.length; ++j) {
				c0 = Math.abs (matrix[i][j]);

				if (c0 > c1) {
					c1 = c0;
				}
			}

			c[i] = c1;
		}

		k = 0;

		for (int j = 0; j < (matrix.length - 1); ++j) {
			pi1 = 0;

			for (int i = j; i < matrix.length; ++i) {
				pi0 = Math.abs (matrix[index[i]][j]);
				pi0 /= c[index[i]];

				if (pi0 > pi1) {
					pi1 = pi0;
					k = i;
				}
			}

			itmp = index[j];
			index[j] = index[k];
			index[k] = itmp;

			for (int i = (j + 1); i < matrix.length; ++i) {
				pj = matrix[index[i]][j] / matrix[index[j]][j];
				matrix[index[i]][j] = pj;

				for (int l = (j + 1); l < matrix.length; ++l) {
					matrix[index[i]][l] -= pj * matrix[index[j]][l];
				}
			}
		}
	}

public static float Determinant(float[][]matrix) {

float temporary[][];
		float result = 0;

		if (matrix.length == 1) {
			result = matrix[0][0];
			return (result);
		}

		if (matrix.length == 2) {
			result = ((matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]));
			return (result);
		}

		for (int i = 0; i < matrix[0].length; i++) {
			temporary = new float[matrix.length - 1][matrix[0].length - 1];

			for (int j = 1; j < matrix.length; j++) {
				for (int k = 0; k < matrix[0].length; k++) {
					if (k < i) {
						temporary[j - 1][k] = matrix[j][k];
					} else if (k > i) {
						temporary[j - 1][k - 1] = matrix[j][k];
					}
				}
			}

			result += matrix[0][i] * Math.pow (-1, (float) i) * Determinant (temporary);
		}
		return (result);


} // determinant function

public static float[][]getMatrix(int row,int col) {

float[][]A = new float[row][col];


for (int i=0;i<row;i++) { 

	for (int j=0;j<col;j++) { 

	s = JOptionPane.showInputDialog("Enter Element No."+(i+1)+"of Row No."+(i+1));
	A[i][j] = Float.parseFloat(s);	

	System.out.print("\t"+A[i][j]);

	} // j
System.out.println();

} // i

return A;

} // getMatrix Function

public static void showMatrix(float[][]M) {

int row = M.length;
int col = M[0].length;

for (int i=0;i<row;i++) {

for (int j=0;j<col;j++) {

System.out.print("\t"+M[i][j]);

//System.out.print("\tRow No."+(i+1)+"\'s Elements ==>\t"+M[i][j]);

}

System.out.println( );

}

} // show Matrix Function

public static float[][]Multiply_Matrices(float[][]a,float[][]b) {

int rowA = a.length;
int colA = a[0].length;

int rowB = b.length;
int colB = b[0].length;

float[][]result = new float[rowA][colB];


if (colA!=rowB) { System.out.println("Rows of A must equals to the Columans of B"); return null; }

for (int i=0;i<rowA;i++) {

	for (int j=0;j<colB;j++) { 

	float sum=0;
	
		for (int k=0;k<colA;k++) {

		sum+= a[i][k] * b[k][j];

		} // k

	result[i][j] = sum;

	} // j

} // i

return result;

} // Multiply Matrices Function

public static float[][]Multiply_Scalar(float[][]a,int scalar) {

int row = a.length;
int col = a[0].length;

float[][]result = new float[row][col];

for (int i=0;i<row;i++) {

for (int j=0;j<col;j++) {

result[i][j] = scalar * a[i][j];

} // j

} // i

return result;

} // multiply with scalar

public static float[][]Transpose_Matrix(float[][]a) {

int row = a.length;
int col = a[0].length;

if (row>col) {

float[][]result = new float[row][col+1];

for (int j=0;j<col+1;j++) {

for (int i=0;i<row-1;i++) {
	result[i][j] = a[j][i];
} // i

} // j

float[][]res = new float[result.length-1][result[0].length];

for (int i=0;i<result.length-1;i++) {

for (int j=0;j<result[0].length;j++) {

//System.out.print("\t"+result[i][j]);

res[i][j] = result[i][j];
} // col
System.out.println();

} // row

 return res;

} // row>col (3x2)

else if (col>row) {

float[][]result = new float[row+1][col];

for (int j=0;j<col-1;j++) {

for (int i=0;i<row+1;i++) {
	result[i][j] = a[j][i];
} // i

} // j

float[][]res = new float[result.length][result[0].length-1];

for (int i=0;i<result.length;i++) {

for (int j=0;j<result[0].length-1;j++) {

//System.out.print("\t"+result[i][j]);

res[i][j] = result[i][j];

} // col
System.out.println();

} // row

return res;

} // col>row (2x3)

else {

float[][]result = new float[row][col];

for (int j=0;j<col;j++) {

for (int i=0;i<row;i++) {

	result[i][j] = a[j][i];

} // i

} // j

return result;

}

} // Transpose Matrix

public static float[][]Adding_Matrices(float[][]a,float[][]b) {

int rowA = a.length;
int colA = a[0].length;

int rowB = b.length;
int colB = b[0].length;

float[][]result = new float[rowA][colA];

if (rowA!=rowB && colA!=colB) { System.out.println("Matrices have to be in same order!!!"); return null; }

for (int i=0;i<rowA;i++) {

for (int j=0;j<colA;j++) {

result[i][j] = a[i][j]+b[i][j];

}

}

return result;

} // adding function

public static float[][]Substracting_Matrices(float[][]a,float[][]b) {

int rowA = a.length;
int colA = a[0].length;

int rowB = b.length;
int colB = b[0].length;

float[][]result = new float[rowA][colA];

if (rowA!=rowB && colA!=colB) { System.out.println("Matrices have to be in same order!!!"); return null; }

for (int i=0;i<rowA;i++) {

for (int j=0;j<colA;j++) {

result[i][j] = a[i][j]-b[i][j];

}

}

return result;

} // substracting function


}