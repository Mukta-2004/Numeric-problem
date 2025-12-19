# Table of Contents

- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
  - [Matrix Inversion](#matrix-inversion)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)

- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)
  - [Secant Method](#secant-method)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)
  - [Newton Raphson Method](#newton-raphson-method)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)

- [Solution of Differential Equations and Differentiation](#solution-of-differential-equations-and-differentiation)
  - [Runge-Kutta Method](#runge-kutta-method)
    - [Theory](#runge-kutta-theory)
    - [Code](#runge-kutta-code)
    - [Input](#runge-kutta-input)
    - [Output](#runge-kutta-output)
  - [Differentiation Method](#differentiation-method)
    - [Theory](#differentiation-theory)
    - [Code](#differentiation-code)
    - [Input](#differentiation-input)
    - [Output](#differentiation-output)

- [Solution of Curve Fitting](#solution-of-curve-fitting)
  - [Linear Equation](#linear-equation-method)
    - [Theory](#linear-equation-theory)
    - [Code](#linear-equation-code)
    - [Input](#linear-equation-input)
    - [Output](#linear-equation-output)
  - [Polynomial Equation](#polynomial-equation-method)
    - [Theory](#polynomial-equation-theory)
    - [Code](#polynomial-equation-code)
    - [Input](#polynomial-equation-input)
    - [Output](#polynomial-equation-output)
  - [Transcendental Equation](#transcendental-equation-method)
    - [Theory](#transcendental-equation-theory)
    - [Code](#transcendental-equation-code)
    - [Input](#transcendental-equation-input)
    - [Output](#transcendental-equation-output)

- [Solution of Interpolation and Approximation](#solution-of-interpolation-and-approximation)
  - [Newton's Forward Interpolation](#newton-forward-interpolation)
    - [Theory](#newton-forward-interpolation-theory)
    - [Code](#newton-forward-interpolation-code)
    - [Input](#newton-forward-interpolation-input)
    - [Output](#newton-forward-interpolation-output)
  - [Newton's Backward Interpolation](#newton-backward-interpolation)
    - [Theory](#newton-backward-interpolation-theory)
    - [Code](#newton-backward-interpolation-code)
    - [Input](#newton-backward-interpolation-input)
    - [Output](#newton-backward-interpolation-output)
  - [Newton's Divided Difference Interpolation](#newton-divided-difference-interpolation)
    - [Theory](#newton-divided-difference-interpolation-theory)
    - [Code](#newton-divided-difference-interpolation-code)
    - [Input](#newton-divided-difference-interpolation-input)
    - [Output](#newton-divided-difference-interpolation-output)

- [Solution of Numerical Integration](#solution-of-numerical-integration)
  - [Simpson's One-Third Rule](#simpson-one-third-rule)
    - [Theory](#simpson-one-third-rule-theory)
    - [Code](#simpson-one-third-rule-code)
    - [Input](#simpson-one-third-rule-input)
    - [Output](#simpson-one-third-rule-output)
  - [Simpson's Three-Eighth Rule](#simpson-three-eighth-rule)
    - [Theory](#simpson-three-eighth-rule-theory)
    - [Code](#simpson-three-eighth-rule-code)
    - [Input](#simpson-three-eighth-rule-input)
    - [Output](#simpson-three-eighth-rule-output)

---


### Solution of Linear Equations

### Gauss Elimination Method

#### Gauss Elimination Theory
[Add your theory content here]

#### Gauss Elimination Code
```python
# Add your code here
```

#### Gauss Elimination Input
```
[Add your input format here]
```

#### Gauss Elimination Output
```
[Add your output format here]
```

---

### Gauss Jordan Elimination Method

#### Gauss Jordan Theory
Introduction:

The Gauss–Jordan Method is a numerical technique used to solve a system of linear equations. It transforms the given system into an equivalent system whose solution can be obtained directly. The method operates on the augmented matrix of the system and reduces it to reduced row echelon form.

Working Principle:

The system of linear equations is first written in augmented matrix form. The method then performs a sequence of elementary row operations to make the matrix diagonal with all diagonal elements equal to one.

For each column:
1. A suitable pivot element is selected.
2. The pivot row is normalized.
3. All other elements in the pivot column are made zero.
4. This process continues until the matrix is reduced. The solution is then read directly from the last column of the matrix.

Special Cases:
1. If a row becomes all zeros except the last column, the system has no solution.
2. If one or more variables do not have a leading one, the system has infinitely many solutions.
3. If every variable has a leading one, the system has a unique solution.

Advantages:
1. Directly provides the solution without back substitution.
2. Systematic and easy to implement.
3. Can detect no solution or infinite solutions.

Limitations:
1. Computationally expensive for large systems.
2. Sensitive to rounding errors.
3. Not efficient compared to iterative methods for very large matrices.

Best and Worst Use Cases:

Works Best When:
1. The number of equations is small to moderate
2. Coefficient matrix is well-conditioned

Works Worst When:
1. The system is very large.
2. The matrix is nearly singular.
3. High precision is required.

Conclusion:

The Gauss–Jordan Method is a powerful direct method for solving systems of linear equations. While simple and reliable, it becomes computationally expensive for large systems and requires careful handling of numerical errors.

#### Gauss Jordan Code
```python
#include<bits/stdc++.h>
#include<fstream>
using namespace std;

const double E = 1e-9;
void print(const vector<vector<double>>& a, ofstream &out)
{
    int n=a.size();
    out<<fixed<<setprecision(3);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<=n;j++)
        {
            out<<setw(8)<<a[i][j]<<" ";
        }
        out<<endl;
    }
    out<<endl;
}

int main()
{
    string inputfile = "InputGaussJordan.txt";
    string outputfile = "OutputGaussJordan.txt";
    ifstream in(inputfile);
    ofstream out(outputfile);

    if(!in)
    {
        cout<<"Input file not found"<<endl;
        return 0;
    }
    if(!out)
    {
        cout<<"Output file not found"<<endl;
        return 0;
    }
    int n;
    in>>n;
    vector<vector<double>> a(n, vector<double>(n + 1));

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<=n;j++)
        {
            in>>a[i][j];
        }
    }
    out<<"Initial Matrix : "<<endl;
    print(a, out);

    for(int col=0,row=0;col<n&&row<n;col++)
    {
        int pivot = row;
        for(int i=row;i<n;i++)
        {
            if(fabs(a[i][col])>fabs(a[pivot][col]))
                pivot=i;
        }
        if (fabs(a[pivot][col]) < E)
            continue;
        swap(a[row],a[pivot]);
        double div = a[row][col];
        for(int j=0;j<=n;j++)
        {
            a[row][j]/=div;
        }

        for(int i=0;i<n;i++)
        {
            if(i!=row)
            {
                double factor = a[i][col];
                for(int j=0;j<=n;j++)
                {
                    a[i][j]-=factor*a[row][j];
                }
            }
        }
        out<<"After processing column "<<col + 1<<endl;
        print(a, out);
        row++;
    }

    int rA = 0, rAug = 0;

    for(int i=0;i<n;i++)
    {
        bool a0 = true;
        bool aug0 = true;
        for(int j=0;j<n;j++)
        {
            if(fabs(a[i][j])>E)
            {
                a0 = false;
                break;
            }
        }
        if(fabs(a[i][n])>E)
            aug0 = false;
        if(!a0)
            rA++;
        if(!a0 || !aug0)
            rAug++;
    }

    if(rA<rAug)
    {
        out<<"Result: No solution"<<endl;
        return 0;
    }
    if(rA<n)
    {
        out<<"Result: Infinitely many solutions"<<endl;
        return 0;
    }
    out<<"Result: Unique solution"<<endl;
    out<<"Solution:"<<endl;
    for(int i=0;i<n;i++)
    {
        out<<"x"<<i + 1<<" = "<<a[i][n]<<endl;
    }
    return 0;
}

```

#### Gauss Jordan Input
```
2
2 1 5
1 -1 1

```

#### Gauss Jordan Output
```
Initial Matrix : 
   2.000    1.000    5.000 
   1.000   -1.000    1.000 

After processing column 1
   1.000    0.500    2.500 
   0.000   -1.500   -1.500 

After processing column 2
   1.000    0.000    2.000 
  -0.000    1.000    1.000 

Result: Unique solution
Solution:
x1 = 2.000
x2 = 1.000

```

---

### LU Decomposition Method

#### LU Decomposition Theory
[Add your theory content here]

#### LU Decomposition Code
```python
# Add your code here
```

#### LU Decomposition Input
```
[Add your input format here]
```

#### LU Decomposition Output
```
[Add your output format here]
```

---

### Matrix Inversion

#### Matrix Inversion Theory
Introduction:

The Matrix Inversion Method is a numerical technique used to find the inverse of a square matrix. Once the inverse of a matrix is obtained, it can be used to solve a system of linear equations of the form A*X = B by computing

X = inverse(A) * B.

Working Principle:

The method forms an augmented matrix by placing the identity matrix beside the given square matrix. Using elementary row operations, the original matrix is transformed into the identity matrix while applying the same operations to the identity matrix.

Steps involved:
1.	Start with the augmented matrix [A | I]
2. Convert matrix A into the identity matrix using row operations
3. The transformed identity matrix becomes the inverse of A
4. If at any stage a pivot element becomes zero, the matrix does not have an inverse.

Special Cases:
1. A matrix with zero determinant is not invertible.
2. If a pivot element is zero, the inversion process fails.
3. Only square matrices can have inverses.

Advantages:
1. Direct method to compute matrix inverse.
2. Useful for solving multiple systems with the same coefficient matrix.

Limitations:
1. Computationally expensive for large matrices.
2. Sensitive to rounding errors.
3. Not suitable for nearly singular matrices.

Best and Worst Use Cases:

Works Best When:
1. Matrix size is small to moderate.
2. Matrix is well-conditioned and non-singular.
Works Worst When:
1.Matrix is very large.
2.Matrix is singular or nearly singular.

Conclusion:

The Matrix Inversion Method is an effective technique for finding the inverse of a matrix and solving linear systems. While straightforward, it becomes inefficient for large matrices and requires careful handling of numerical precision.

#### Matrix Inversion Code
```python
#include<bits/stdc++.h>
#include<fstream>
using namespace std;

bool inmat(vector<vector<double>>& A, vector<vector<double>>& inv)
{
    int n = A.size();
    inv.assign(n,vector<double>(n, 0));
    for(int i=0; i<n; i++)
    {
        inv[i][i]=1;
    }
    for(int i=0; i<n; i++)
    {
        double pivot=A[i][i];
        if(pivot==0)
        {
            return false;
        }
        for(int j=0;j<n;j++)
        {
            A[i][j]/=pivot;
            inv[i][j]/=pivot;
        }
        for(int k=0;k<n;k++)
        {
            if(k!=i)
            {
                double factor=A[k][i];
                for(int j=0;j<n;j++)
                {
                    A[k][j]-=factor*A[i][j];
                    inv[k][j]-=factor*inv[i][j];
                }
            }
        }
    }

    return true;
}

int main()
{
    string inputfile = "InputMatrixInverse.txt";
    string outputfile = "OutputMatrixInverse.txt";
    ifstream in(inputfile);
    if (!in)
    {
        cout << "Input file error" << endl;
        return 1;
    }

    int n;
    in>>n;
    vector<vector<double>> mat(n,vector<double>(n));
    for (int i=0;i<n;i++)
    {
        for (int j = 0; j < n; j++)
        {
            in>>mat[i][j];
        }
    }
    in.close();

    ofstream out(outputfile);
    if (!out)
    {
        cout<<"Output file error"<<endl;
        return 1;
    }
    vector<vector<double>> inv;

    if (!inmat(mat, inv))
    {
        out<<"Matrix is not invertible"<<endl;
    }
    else
    {
        out<<"Inverse Matrix : "<<endl;
        out<<fixed<<setprecision(6);
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                out<<inv[i][j]<<" ";
            }
            out<<endl;
        }
    }
    out.close();
    return 0;
}


```

#### Matrix Inversion Input
```
3
2 1 1
1 3 2
1 0 0

```

#### Matrix Inversion Output
```
Inverse Matrix : 
0.000000 0.000000 1.000000 
-2.000000 1.000000 3.000000 
3.000000 -1.000000 -5.000000 

```

---

### Solution of Non-Linear Equations

### Bisection Method

#### Bisection Theory
[Add your theory content here]

#### Bisection Code
```python
# Add your code here
```

#### Bisection Input
```
[Add your input format here]
```

#### Bisection Output
```
[Add your output format here]
```

---

### False Position Method

#### False Position Theory
[Add your theory content here]

#### False Position Code
```python
# Add your code here
```

#### False Position Input
```
[Add your input format here]
```

#### False Position Output
```
[Add your output format here]
```
### Secant Method
#### Secant Theory
[Add your theory content here]
#### Secant Code
```python
# Add your code here
```
#### Secant Input
```
[Add your input format here]
```
#### Secant Output
```
[Add your output format here]
```

### Newton Raphson Method
#### Newton Raphson Theory
[Add your theory content here]
#### Newton Raphson Code
```python
# Add your code here
```
#### Newton Raphson Input
```
[Add your input format here]
```
#### Newton Raphson Output
```
[Add your output format here]
```

### Solution of Differential Equations and Differentiation
### Runge-Kutta Method
#### Runge-Kutta Theory
[Add your theory content here]
#### Runge-Kutta Code
```python
# Add your code here
```
#### Runge-Kutta Input
```
[Add your input format here]
```
#### Runge-Kutta Output
```
[Add your output format here]
```

### Differentiation Method
#### Differentiation Theory
[Add your theory content here]
#### Differentiation Code
```python
# Add your code here
```
#### Differentiation Input
```
[Add your input format here]
```
#### Differentiation Output
```
[Add your output format here]
```

---

## Solution of Curve Fitting

### Linear Equation

#### Linear Equation Theory
Introduction:

Linear Regression is a numerical method used to determine the best fitting straight line between two variables x and y. The relationship is represented by a linear equation:

y = a + b*x

where a is the intercept and b is the slope of the line.

Working Principle:

The method uses the Least Squares approach to minimize the total squared error between observed data points and estimated values. Given n data points, required sums are calculated from x and y values. The slope and intercept are obtained using:

b = (n*sum(x*y) - sum(x)sum(y)) / (n*sum(x*x) - (sum(x))^2) ; 

a = (sum(y) - b*sum(x)) / n

Using these values, the best fit linear equation is formed.

Special Cases:
1. If all x values are equal, the denominator becomes zero and regression is not possible.
2. At least two data points are required.
3. If data points lie exactly on a straight line, the regression gives a perfect fit.

Advantages:
1. Simple and easy to implement.
2. Requires less computation.
3. Works well for large datasets.
4. Useful for prediction and trend analysis.

Best and Worst Use Cases:

Works Best When:
1. Data follows an approximately linear pattern
2. Outliers are minimal
   
Works Worst When:
1. Data is highly nonlinear
2. Extreme outliers are present
3. Very small datasets

Conclusion:
Linear Regression is a fundamental numerical technique for modeling linear relationships. While easy to apply and efficient, its accuracy depends on the nature and quality of the data.

#### Linear Equation Code
```python
#include<bits/stdc++.h>
#include <fstream>
using namespace std;

int main()
{
    string ifile = "InputLinear.txt";
    string ofile = "OutputLinear.txt";

    ifstream in(ifile);
    if (!in)
    {
        cout << "Error opening input file"<<endl;
        return 1;
    }
    int n;
    in >> n;
    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++)
    {
        in >> x[i];
    }
    for (int i = 0; i < n; i++)
    {
        in >> y[i];
    }
    in.close();
    double sumx = 0;
    double sumy = 0;
    double sumxy = 0;
    double sumx2 = 0;

    for (int i = 0; i < n; i++)
    {
        sumx+= x[i];
        sumy+= y[i];
        sumxy+= x[i]*y[i];
        sumx2+= x[i]*x[i];
    }
    double b = (n*sumxy-sumx*sumy)/(n*sumx2-sumx*sumx);
    double a = (sumy-b*sumx)/n;

    ofstream out(ofile);
    if (!out)
    {
        cout << "Error opening output file"<<endl;
        return 1;
    }
    out << fixed << setprecision(4);
    out << "Number of data points: " << n <<endl;
    out << "x values:\n";
    for (int i = 0; i < n; i++)
    {
        out << x[i] << " ";
    }
    out <<endl;
    out << "y values:\n";
    for (int i = 0; i < n; i++)
    {
        out << y[i] << " ";
    }
    out << endl;
    out << "y = " << a << " + " << b << "x"<<endl;
    out.close();
    return 0;
}

```

#### Linear Equation Input
```
5
1 2 3 4 5
2 4 5 4 5

```

#### Linear Equation Output
```
Number of data points: 5
x values:
1.0000 2.0000 3.0000 4.0000 5.0000 
y values:
2.0000 4.0000 5.0000 4.0000 5.0000 
y = 2.2000 + 0.6000x

```

---

### Polynomial Equation

#### Polynomial Equation Theory
[Add your theory content here]

#### Polynomial Equation Code
```python
# Add your code here
```

#### Polynomial Equation Input
```
[Add your input format here]
```

#### Polynomial Equation Output
```
[Add your output format here]
```

---

### Transcendental Equation
#### Transcendental Equation Theory
Introduction:

Transcendental Regression is used when the relationship between variables is non-linear but can be converted into a linear form using mathematical transformations. This method applies linear regression after transforming the given data.
In this program, two models are supported:

Power model: y = a*x^b

Exponential model: p = p0 * e^(k*t)


Working Principle:

The given nonlinear equations are converted into linear equations using logarithms.

For the power model: log(y) = log(a) + b * log(x)

For the exponential model: log(p) = log(p0) + k * t

After transformation, the variables become linear and standard least squares linear regression is applied. The slope and intercept are calculated, and the original constants are recovered using exponential functions.

Special Cases:
1. All y values must be positive because logarithm of zero or negative numbers is undefined.
2. For the power model, x values must also be positive.
3. At least two data points are required.
4. If transformed data is perfectly linear, the regression gives an exact fit.

Advantages:
1. Allows regression for nonlinear models.
2. Uses simple linear regression after transformation.

Best and Worst Use Cases:

Works Best When:

1. Data follows power or exponential behavior.
2. All data values are positive.
3. Noise in data is minimal.
   
Works Worst When:

1. Data contains zero or negative values
2. Relationship cannot be linearized
3. Dataset is very small or highly noisy

Conclusion:

Transcendental Regression extends linear regression to nonlinear equations by using transformations. It is effective for power and exponential models but requires careful handling of data values and assumptions.
#### Transcendental Equation Code
```python
#include<bits/stdc++.h>
#include<fstream>
using namespace std;

int main()
{
    string inputFile = "InputTranscendental.txt";
    string outputFile = "OutputTranscendental.txt";

    ifstream in("InputTranscendental.txt");
    if (!in )
    {
        cout << "Input file error"<<endl;
        return 1;
    }

    int choice, n;
    in >> choice;// 1-> y = ax^b ; 2-> p = p0*e^(kt)
    in >> n;
    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++)
    {
        in >> x[i];
    }
    for (int i = 0; i < n; i++)
    {
        in >> y[i];
    }
    in.close();

    double sx = 0;
    double sy = 0;
    double sxx = 0;
    double sxy = 0;

    for (int i = 0; i < n; i++)
    {
        double X, Y;
        Y = log(y[i]);
        if (choice == 1)
        {
            X = log(x[i]);
        }
        else if(choice == 2)
        {
            X = x[i];
        }
        else
        {
            cout << "Invalid choice"<<endl;
        }
        sx  += X;
        sy  += Y;
        sxx += X*X;
        sxy += X*Y;
    }

    double b = (n*sxy-sx*sy)/(n*sxx-sx*sx);
    double A = (sy-b*sx)/n;
    ofstream out(outputFile);
    if (!out)
    {
        cout << "Output file error"<<endl;
        return 1;
    }
    out << fixed << setprecision(6);
    out << "Number of data points: " << n <<endl;
    out << "x values:"<<endl;
    for (int i = 0; i < n; i++)
    {
        out << x[i] << " ";
    }
    out <<endl;
    out << "y values:"<<endl;
    for (int i = 0; i < n; i++)
    {
        out << y[i] << " ";
    }
    out<<endl;
    if (choice == 1)
    {
        double a = exp(A);
        out << "y = " << a << " * x^(" << b << ")"<<endl;
    }
    else if (choice == 2)
    {
        double p0 = exp(A);
        double k  = b;
        out << "p = " << p0 << " * e^(" << k << " t)"<<endl;
    }
    else
    {
        out << "Invalid choice"<<endl;
    }
    out.close();
    return 0;
}


```
#### Transcendental Equation Input
```
1
4
1 2 3 4
2.7 7.4 20.1 54.6

```
#### Transcendental Equation Output
```
Number of data points: 4
x values:
1.000000 2.000000 3.000000 4.000000 
y values:
2.700000 7.400000 20.100000 54.600000 
y = 2.276154 * x^(2.109951)

```
### Solution of Interpolation and Approximation
### Newton's Forward Interpolation
#### Newton's Forward Interpolation Theory
[Add your theory content here]
#### Newton's Forward Interpolation Code
```python
#Add your code here
```
#### Newton's Forward Interpolation Input
```
[Add your input format here]
```
#### Newton's Forward Interpolation Output
```
[Add your output format here]
```

### Newton's Backward Interpolation
#### Newton's Backward Interpolation Theory
[Add your theory content here]
#### Newton's Backward Interpolation Code
```python
#Add your code here
```
#### Newton's Backward Interpolation Input
```
[Add your input format here]
```
#### Newton's Backward Interpolation Output
```
[Add your output format here]
```

### Newton's Divided Difference Interpolation
#### Newton's Divided Difference Interpolation Theory
[Add your theory content here]
#### Newton's Divided Difference Interpolation Code
```python
#Add your code here
```
#### Newton's Divided Difference Interpolation Input
```
[Add your input format here]
```
#### Newton's Divided Difference Interpolation Output
```
[Add your output format here]
```

### Solution of Numerical Integration
### Simpson's One-Third Rule
#### Simpson's One-Third Rule Theory
[Add your theory content here]
#### Simpson's One-Third Rule Code
```python
#Add your code here
```
#### Simpson's One-Third Rule Input
```
[Add your input format here]
```
#### Simpson's One-Third Rule Output
```
[Add your output format here]
```
### Simpson's Three-Eighth Rule
#### Simpson's Three-Eighth Rule Theory
[Add your theory content here]
#### Simpson's Three-Eighth Rule Code
```python
#Add your code here
```
#### Simpson's Three-Eighth Rule Input
```
[Add your input format here]
```
#### Simpson's Three-Eighth Rule Output
```
[Add your output format here]
```



---





