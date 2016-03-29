#include <iostream>
#include <math.h>

using namespace std;

double Cosinus(double, double);
double Sinus(double, double);
double Exponent(double, double);
double** AllocMemory(int);
void EllementsInMyMatrix(double**, int, double);
void EllementsInMatrix(double**, int);
void DifferenceMatrix(double**, double**, double**, int);
void Display(double**, int);
void FreeMemory(double**, int);
double MaxElementInDifference(double**, int);



int main()
{
	int n;
	double eps;
	while (true)
	{
		cout << "Please, enter size of matrix: ";
		cin >> n;
		cout << "Please, enter epsilon: ";
		cin >> eps;
		/*system("cls");*/
		if ((eps > 0 && eps < 1) || (n > 1)) break;
		cout << "Error! Please, enter other size or value of epsilon " << endl;
	}
	double** a = AllocMemory(n);
	EllementsInMyMatrix(a, n, eps);
	cout << "Matrix with our function: " << endl;
	Display(a, n);
	double** b = AllocMemory(n);
	cout << endl << "Matrix with standart function: " << endl;
	EllementsInMatrix(b, n);
	Display(b, n);
	double** c = AllocMemory(n);
	DifferenceMatrix(a, b, c, n);
	cout << endl << "Difference:" << endl;
	Display(c, n);
	double max = MaxElementInDifference(c, n);
	cout << endl << "Max ellement: " << endl;
	cout << max << endl;
	system("pause");
	FreeMemory(a, n);
	FreeMemory(b, n);
	FreeMemory(c, n);
	return 0;
}

double Cosinus(double a, double eps)
{
	int m = 1;
	double sum = 1, n = (-1 * a * a) / (m * (m + 1));
	while (fabs(n) >= eps)
	{
		m += 2;
		sum += n;
		n *= (-1 * a * a) / (m * (m + 1));
	}
	return sum;
}

double Sinus(double x, double eps)
{
	double n = x;
	double sum = 0.0;
	int i = 1;

	while (fabs(n) > eps)
	{
		sum += n;
		n *= -1.0 * x * x / ((2 * i) * (2 * i + 1));
		i++;
	}
	return sum;
}

double Exponent(double x, double eps)
{
	double sum = 0, p = 1;
	int i = 1;
	while (fabs(p) > eps)
	{
		sum += p;
		p *= Sinus(x, eps) / i;
		i++;
	}
	return sum;
}

double** AllocMemory(int n)
{
	double** a = new double*[n];
	for (int i = 0; i < n; i++)
		a[i] = new double[n];
	return a;

}

void EllementsInMyMatrix(double** a, int n, double eps)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a[i][j] = ((Exponent((i + j), eps) + (Cosinus((i + j),eps)*Cosinus((i + j), eps))) / (Sinus((i + 1), eps)*Sinus((i + 1), eps)));
			if (i == j)
				a[i][j] = 0;

		}
	}
}

void EllementsInMatrix(double** a, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a[i][j] = (exp(sin(i + j)) + (cos(i + j)*cos(i + j))) / (sin(i + 1)*sin(i + 1));
			if (i == j)
				a[i][j] = 0;

		}
	}
}

void DifferenceMatrix(double** a, double** b, double** c, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			c[i][j] = fabs(a[i][j] - b[i][j]);
		}
	}
}

void Display(double** a, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout.width(15);
			cout << a[i][j];
		}
		cout << endl;
	}
}

void FreeMemory(double** a, int n)
{
	for (int i = 0; i < n; i++)
		delete[] a[i];
	delete[] a;
}

double MaxElementInDifference(double** c, int n)
{
	double max = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (c[i][j] > max) max = c[i][j];
		}
	}
	return max;
}
