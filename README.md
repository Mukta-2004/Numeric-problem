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
Gauss Elimination Method

The Gauss elimination method is a fundamental numerical technique used to solve systems of linear equations. It converts a system into an upper triangular form and then solves it using back substitution, simplifying the solution process. This method is applicable to square and rectangular systems and forms the basis for many advanced numerical methods in linear algebra.

Problem Definition
We aim to solve a system of linear equations in matrix form:

𝐴𝑥=𝑏

where A is the coefficient matrix, x is the vector of unknowns, and b is the constants vector.
The method proceeds in two main phases:

1. Forward Elimination Phase:

In this phase, the augmented matrix [A∣b] is transformed into an upper triangular matrix. This is done by eliminating the entries below the diagonal in each column using the following steps:
Select a pivot element in the current column.
Use the pivot to eliminate all entries below it, replacing each row with a combination of itself and the pivot row.
Repeat the process for all columns until the augmented matrix has zeros below the diagonal.
The result is an upper triangular matrix, where the system of equations can be solved more easily.

3. Back Substitution Phase:
   
Once the system is in upper triangular form, the values of the unknowns are determined starting from the last row upwards:
Solve the last equation for the last unknown.
Substitute this value into the above equations to solve for other unknowns iteratively.
Continue this process until all unknowns are computed.
This phase uses simple algebraic substitution and is straightforward once the matrix is triangular.

Applications
1. Solving linear algebraic systems in engineering and physics.
2. Finding solutions for structural analysis, electrical circuits, and mechanics problems.
3. Forms the foundation for more advanced methods like LU decomposition and matrix inversion.


Advantages:

1. Conceptually simple and easy to implement.
2. Systematic approach that works for most non-singular systems.
3. Forms a foundation for more advanced numerical methods.

Disadvantages:

1. Computationally expensive for large systems (O(n³) complexity).
2. Can suffer from round-off errors in floating-point arithmetic.
3. Requires row swaps (partial pivoting) to avoid division by zero or improve numerical stability.




#### Gauss Elimination Code
```cpp
#include <bits/stdc++.h>
#include<fstream>
using namespace std;

int main(){
  string inputFile = "Gauss Elimination input.txt";
    string outputFile ="Gauss Elimination output.txt";
    ifstream in(inputFile);

if(!in){
    cout<<" Input file error"<<endl;
    return 1;
}

int n;
in >> n;
vector<vector<double>> a(n, vector<double>(n + 1));
vector<double> x(n);
 for (int i = 0; i < n; i++){
    for (int j = 0; j <= n; j++){
            in >> a[i][j];
    }
 }
 in.close();
ofstream out(outputFile);
 if(!out){
    cout<<" Output file error"<<endl;
    return 1;
}
out << "Number of unknowns: "<<n << endl;
for (int i = 0; i < n; i++){
  if (a[i][i] == 0){
    out << "Mathematical Error" << endl;
    return 0;
  }
  for (int j = i + 1; j < n; j++){
    double rtio = (a[j][i] / a[i][i]);
    for (int k = 0; k <= n; k++){
         a[j][k] = a[j][k] - (rtio * a[i][k]);
    }
  }
}

 for (int i = n - 1; i >= 0; i--){
    x[i] = a[i][n];
    for (int j = i + 1; j < n; j++){
            x[i] = x[i] - a[i][j] * x[j];
    }
        x[i] = x[i] / a[i][i];
 }

out << "Solution:" << endl;
for (int i = 0; i < n; i++){
        out << "x" << i + 1 << " = "<< fixed << setprecision(6) << x[i] << endl;
}

out.close();
return 0;
}

```

#### Gauss Elimination Input
```
2
3 2 16
4 -1 9

```

#### Gauss Elimination Output
```
Number of unknowns: 2
Solution:
x1 = 3.090909
x2 = 3.363636

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
```cpp
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
```cpp
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
Introduction:

The False Position Method is a numerical technique used to find the root of a nonlinear equation f(x) = 0. It is a bracketing method, meaning it starts with an interval where the function values at the endpoints have opposite signs, ensuring that a root lies between them.

Working Principle:

The method begins with two initial guesses a and b such that f(a) and f(b) have opposite signs. A new approximation c is calculated by drawing a straight line between the points (a, f(a)) and (b, f(b)) and finding where it intersects the x-axis. The approximation is computed using:

c = b - f(b)*(b - a) / (f(b) - f(a))

The value of c replaces either a or b depending on the sign of f(c), so the root remains bracketed. This process is repeated until the function value at c is smaller than the given tolerance or the maximum number of iterations is reached.

Special Cases:
1. The method fails if f(a) and f(b) have the same sign.
2. If one endpoint remains fixed, convergence may become slow.
3. The method gives an exact root if f(c) becomes zero.

Advantages:
1. Always maintains the root within the interval.
2. Guaranteed convergence if initial condition is satisfied.
3. More accurate than the bisection method in many cases.

Limitations:
1. Convergence can be slow for certain functions.
2. One endpoint may remain unchanged for many iterations.
3. Slower compared to open methods like Newton-Raphson.

Best and Worst Use Cases:

Works Best When:
1. f(a) and f(b) have opposite signs.
2. The function is continuous in the interval.
Works Worst When:
1. Function slope is very small near one endpoint.
2. Poor initial interval selection.

Conclusion:

The False Position Method is a reliable root-finding technique that combines the safety of bracketing methods with faster convergence than bisection. However, its performance depends on the nature of the function and interval selection.

#### False Position Code
```cpp
#include<bits/stdc++.h>
#include <fstream>
using namespace std;

double f(double x,double a4,double a3,double a2,double a1,double a0)
{
    return a4*pow(x,4)+a3*pow(x,3)+a2*pow(x,2)+a1*x+a0;
}
double fp(double a,double b,double E,int mIter,double a4,double a3,double a2,double a1,double a0,ofstream &out)
{
    if (f(a,a4,a3,a2,a1,a0)*f(b,a4,a3,a2,a1,a0)>=0)
    {
        out<<"Invalid as f(a) and f(b) must have opposite signs"<<endl;;
        return 1;
    }

    double c;
    for (int i = 1; i <= mIter; i++)
    {
        c = a-(f(a,a4,a3,a2,a1,a0)*(b-a))/(f(b,a4,a3,a2,a1,a0)-f(a,a4,a3,a2,a1,a0));
        out<<"Iteration "<<i<<" , c = "<<c<<" , f(c) = "<<f(c,a4,a3,a2,a1,a0)<<endl;

        if(fabs(f(c,a4,a3,a2,a1,a0))<E)
        {
            out<<"Converged after "<<i<<" iterations"<<endl;
            return c;
        }
        if(f(a,a4,a3,a2,a1,a0)*f(c,a4,a3,a2,a1,a0)<0)
            b = c;
        else
            a = c;
    }

    out<<"Max iterations reached"<<endl;
    return c;
}

int main()
{
    string inputFile = "InputFP.txt";
    string outputFile = "OutputFP.txt";

    ifstream in(inputFile);
    if (!in)
    {
        cout << "Input file error"<<endl;
        return 1;
    }

    double a4, a3, a2, a1, a0;
    double a, b, E;
    int mIter;

    in>>a4>>a3>>a2>>a1>>a0;
    in>>a>>b;
    in>>E;
    in>>mIter;
    in.close();

    ofstream out(outputFile);
    if (!out)
    {
        cout<<"Output file error"<<endl;
        return 1;
    }
    out<<"Polynomial: "<<a4<<"x^4 + "<<a3<<"x^3 + "<<a2<<"x^2 + "<<a1<<"x + "<<a0<<endl;
    out<<"Interval: ["<<a<<", "<<b<<"]"<<endl;
    out<<"Tolerance: "<<E<< endl;
    out<<"Max Iterations: "<<mIter<<endl;
    double root = fp(a,b,E,mIter,a4,a3,a2,a1,a0,out);

    out<<"Approximate Root = "<<root<<endl;

    out.close();

    return 0;
}


```

#### False Position Input
```
1 -3 0 2 -1
-1 0
0.0001
20

```

#### False Position Output
```
Polynomial: 1x^4 + -3x^3 + 0x^2 + 2x + -1
Interval: [-1, 0]
Tolerance: 0.0001
Max Iterations: 20
Iteration 1 , c = -0.5 , f(c) = -1.5625
Iteration 2 , c = -0.804878 , f(c) = -0.625805
Iteration 3 , c = -0.879984 , f(c) = -0.116009
Iteration 4 , c = -0.89246 , f(c) = -0.0180395
Iteration 5 , c = -0.894366 , f(c) = -0.00272592
Iteration 6 , c = -0.894653 , f(c) = -0.000410117
Iteration 7 , c = -0.894696 , f(c) = -6.16619e-05
Converged after 7 iterations
Approximate Root = -0.894696

```
### Secant Method
#### Secant Theory
Introduction:

The Secant Method is a numerical technique used to find the root of a nonlinear equation f(x) = 0. It is an open method, meaning it does not require the root to be bracketed between two values. Instead, it uses two initial approximations to generate successive better estimates.

Working Principle:

The method starts with two initial guesses x0 and x1. A straight line (secant) is drawn between the points (x0, f(x0)) and (x1, f(x1)). The point where this line intersects the x-axis gives the next approximation.
The iteration formula is:

x2 = x1 - f(x1)*(x1 - x0)/(f(x1) - f(x0))

After each iteration, the older value is replaced and the process continues until the function value at the new point is smaller than the given tolerance or the maximum number of iterations is reached.

Special Cases:
1. The method fails if f(x1) - f(x0) becomes zero (division by zero).
2. Poor initial guesses may lead to slow convergence or divergence.
3. Convergence is not guaranteed for all functions.

Advantages:
1. Faster convergence than bisection and false position methods.
2. Does not require derivative calculation.

Limitations:
1. No guarantee of convergence.
2. Sensitive to initial guesses.
3. May diverge for certain functions.

Best and Worst Use Cases:

Works Best When:
1. Initial guesses are close to the actual root.
2. The function is smooth and continuous.
Works Worst When:
1. Initial guesses are far from the root.
2. Function behavior is highly irregular.
3. f(x1) is very close to f(x0).

Conclusion:

The Secant Method is an efficient root-finding technique that improves upon bracketing methods in speed. However, its success strongly depends on the choice of initial guesses.
#### Secant Code
```cpp
#include <bits/stdc++.h>
#include <fstream>
using namespace std;

double f(double x, double a4, double a3, double a2, double a1, double a0)
{
    return a4*pow(x,4) + a3*pow(x,3) + a2*pow(x,2) + a1*x + a0;
}

double secant(double x0,double x1,double E,int mIter,double a4,double a3,double a2,double a1,double a0,ofstream &out)
{
    double x2;
    for(int i = 1; i <= mIter; i++)
    {
        double f0 = f(x0, a4, a3, a2, a1, a0);
        double f1 = f(x1, a4, a3, a2, a1, a0);

        if(fabs(f1 - f0) < 1e-12)
        {
            out<<"Division by zero occured"<<endl;
            return x1;
        }
        x2 = x1-f1*(x1-x0)/(f1-f0);
        out<<"Iteration "<<i<<" , x = "<<x2<<" , f(x) = "<<f(x2, a4, a3, a2, a1, a0)<<endl;

        if(fabs(f(x2, a4, a3, a2, a1, a0))< E)
        {
            out<<"Converged after "<<i<<" iterations"<<endl;
            return x2;
        }
        x0 = x1;
        x1 = x2;
    }
    out<<"Max iterations reached"<<endl;
    return x2;
}

int main()
{
    string inputFile = "InputSecant.txt";
    string outputFile = "OutputSecant.txt";

    ifstream in(inputFile);
    if(!in)
    {
        cout<<"Input file error"<<endl;
        return 1;
    }

    double a4, a3, a2, a1, a0;
    double x0, x1, E;
    int mIter;

    in >> a4 >> a3 >> a2 >> a1 >> a0;
    in >> x0 >> x1;
    in >> E;
    in >> mIter;
    in.close();

    ofstream out(outputFile);
    if(!out)
    {
        cout << "Output file error" << endl;
        return 1;
    }

    out<<"Polynomial: "<<a4<<"x^4 + "<<a3<<"x^3 + "<<a2<<"x^2 + "<<a1<<"x + "<<a0<<endl;
    out<<"Initial guesses: x0 = "<<x0<<", x1 = "<<x1<<endl;
    out<<"Tolerance: "<<E<<endl;
    out<<"Max Iterations: "<<mIter<<endl;

    double root = secant(x0,x1,E,mIter,a4,a3,a2,a1,a0,out);
    out<<"Approximate Root = "<<root<<endl;

    out.close();
    return 0;
}

```
#### Secant Input
```
1 -3 0 2 -1
-1 -0.8
0.0001
20

```
#### Secant Output
```
Polynomial: 1x^4 + -3x^3 + 0x^2 + 2x + -1
Initial guesses: x0 = -1, x1 = -0.8
Tolerance: 0.0001
Max Iterations: 20
Iteration 1 , x = -0.87911 , f(x) = -0.122726
Iteration 2 , x = -0.897371 , f(x) = 0.0216173
Iteration 3 , x = -0.894636 , f(x) = -0.00054223
Iteration 4 , x = -0.894703 , f(x) = -2.29816e-06
Converged after 4 iterations
Approximate Root = -0.894703

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
Numerical differentiation is a technique used to approximate the derivative of a function using discrete data points. In many real-life and engineering problems, the exact mathematical expression of a function may not be available. Instead, the function is known only through experimental observations or tabulated values at equally spaced intervals. In such situations, numerical methods are essential for estimating derivatives.

Newton’s interpolation formulas provide an efficient way to approximate derivatives from tabulated data. By constructing an interpolation polynomial using known data points and then differentiating it, numerical values of derivatives can be obtained. Depending on the location of the point at which the derivative is required, either Newton’s Forward Interpolation Formula or Newton’s Backward Interpolation Formula is applied.

Newton’s Forward Formula is used when the required derivative lies near the beginning of the data table.

Newton’s Backward Formula is used when the derivative is required near the end of the data table.

These methods rely on finite difference tables and are widely used due to their simplicity and effectiveness.

Newton’s Forward Difference Formula for Derivative:

When the derivative is required near the beginning of the table:

dy/dx=1/h[Δy0− 1/2Δ2y0+ 1/3Δ3y0−⋯]

 Newton’s Backward Difference Formula for Derivative:

When the derivative is required near the end of the table:

dy/dx= 1/h[∇yn +1/2∇2yn+ 1/3∇3yn+⋯ ]

Procedure (Step-by-Step)
1.	Collect Data Obtain function values y0,y1,…,yn at equally spaced points x0,x1,…,xn.
Calculate the step size:
ℎ= x1-x0
2. Construct Difference Table
3. Compute u
4. Apply Newton’s Formula for First Derivative
5. Compute the Derivative. Substitute values from the difference table into the formula. Perform the arithmetic carefully to get the approximate derivative at the required point.

Error:
   Relative error is a measure of how large the error is in comparison to the true value. It tells us how significant the error is relative to the size of the quantity being measured.

-> It is dimensionless, meaning it does not depend on the units of measurement.
-> Often expressed as a fraction or percentage.
 
     Relative error = |True value - Approximate value| / True value .

 Applications:

1.Engineering Applications:
Used to calculate velocity, acceleration, and rate of change of physical quantities from measured data.

2.Physics and Mechanics:
Applied in motion analysis, heat transfer, and fluid flow problems where data is experimental.

3.Experimental and Scientific Data Analysis:
Useful when data is obtained from experiments and an analytical function is unavailable.

4.Economics and Statistics:
Used to estimate growth rates, marginal cost, and trend analysis from tabulated data.

5.Computer Science and Numerical Simulations:
Applied in numerical modeling and simulations requiring derivative approximations.


#### Differentiation Code
```cpp
#include <iostream>
#include<bits/stdc++.h>
#include<fstream>

using namespace std;

void BackWard(int n,vector<double>&x, vector<double>&y, double X,double result ,vector<vector<double>>&F,ifstream &in, ofstream &out )
{
    //difference table calculation

    for(int j=1;j<n;j++)
    {
        for(int i=n-1;i>=j;i--)
            F[i][j]=F[i][j-1]-F[i-1][j-1];
    }

out<<"Difference table:\n";
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
            out<<F[i][j]<<"   ";
        out<<endl;
    }

    double h= x[1]-x[0];
    double v=(X-x[n-1])/h;

result= (1/h)*(F[n-1][1]+ ((2*v+1)/2)*F[n-1][2] + ((3*v*v + 6*v +2)/6)*F[n-1][3] + (4*(pow(v,3)+ 18*pow(v,2) + 22*v+6)/24)*F[n-1][4] + (5*(pow(v,4) + 40*pow(v,3) + 105*pow(v,2) +100*pow(v,1) + 24)/120));



    out<<"Backward_result: "<<result<<endl;
}


void ForWard(int n,vector<double>&x, vector<double>&y, double X,double result, vector<vector<double>>&F,ifstream &in, ofstream &out )
{


    //difference table calculation

    for(int j=1;j<n;j++)
    {
        for(int i=0;i<n-j;i++)
            F[i][j]=F[i+1][j-1]-F[i][j-1];
    }
out<<"Difference table:\n";
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
            out<<F[i][j]<<"   ";
        out<<endl;
    }

    double h= x[1]-x[0];
    double u=(X-x[0])/h;

   /* for(int i=0;i<n;i++)
    {
        double key=1;
        for(int j=0;j<i;j++)
        {
            key*=(u-j)/(j+1);
        }

        result+=key*F[0][i];
    }*/
result= (1/h)*(F[0][1]+ ((2*u -1)/2)*F[0][2] + ((3*u*u - 6*u +2)/6)*F[0][3] + (4*(pow(u,3)- 18*pow(u,2) + 22*pow(u,1)- 6)/24)*F[0][4] + (5*(pow(u,4) - 40*pow(u,3) + 105*pow(u,2) -100*pow(u,1) + 24)/120));


    out<<"Forward_result: "<<result<<endl;
}
int main()
{
    string inputFile="Diff_input.txt";
    string outputFile="Diff_output.txt";

    ifstream in(inputFile);
    if(!in)
    {
        cout<<"Input file error!"<<endl;
        return 1;
    }

    ofstream out(outputFile);
    if(!out)
    {
         cout<<"Output file error!"<<endl;
        return 1;
    }

    int n;
    in>>n;
    vector<double>x(n);
    vector<double>y(n);

    for(int i=0;i<n;i++)
        in>>x[i];

    double diff=x[1]-x[0];

    for(int i=2;i<n;i++)
    {
        if((x[i]-x[i-1])!= diff)
           {
               out<<"Difference are not equal\n";
               return 0;
           }
    }

    for(int i=0;i<n;i++)
        in>>y[i];


    vector<vector<double>>F(n,vector<double>(n,0));

    for(int i=0;i<n;i++)
        F[i][0]=y[i];

         double result=0;
    double X;
    in>>X;

    if(X<(x[0]+x[n-1])/2)
        ForWard(n,x,y,X,result,F,in,out);
    else
        BackWard(n,x,y,X,result,F,in,out);

in.close();
out.close();
    return 0;
}

```
#### Differentiation Input
```
4
1 3 5 7
3 5 7 9
8

```
#### Differentiation Output
```
Difference table:
3   0   0   0   
5   2   0   0   
7   2   0   0   
9   2   0   0   
Backward_result: 3.19401

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
```cpp
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
Polynomial curve fitting regression is a statistical technique used to model non-linear relationships between a dependent variable y and an independent variable x. Unlike linear regression, which assumes a constant rate of change, polynomial regression allows the rate of change to vary with 𝑥 by including higher powers of x.

The degree of the polynomial determines the complexity of the curve. A first-degree polynomial represents a straight line, a second-degree polynomial represents a parabolic curve, and higher-degree polynomials can model more complex trends. The primary objective of polynomial regression is to find a polynomial equation that best approximates the given data points by minimizing the overall error.

Polynomial regression is widely used in engineering, physics, economics, and data modeling for trend analysis, prediction, and interpolation where linear models are insufficient.

Procedure of Polynomial Regression:

Step 1: Data Collection: 

Collect the experimental or observed data points: (x1,y1),(x2,y2),…,(xm,ym). These data points represent the independent variable x and the corresponding dependent variable y.

Step 2: Selection of Polynomial Degree :

Choose the degree 𝑛 of the polynomial based on the behavior of the data. The assumed regression model is: y=a0+a1x1+a2x2+⋯+anxn
The choice of n is crucial, as a polynomial of too low a degree may fail to capture the trend, while a polynomial of too high a degree may lead to overfitting.

Step 3: Application of Least Squares Method

To determine the coefficients a0,a1,…,an, the least squares method is employed. This method minimizes the sum of the squares of the residuals (errors) between the observed values 𝑦𝑖 and the predicted values y^i.

The error function is defined as: S=i=1∑m(yi−y^i)2

Step 4: Formation of Normal Equations

Differentiate the error function 𝑆 with respect to each coefficient a0,a1,…,an and equate the derivatives to zero. This results in a system of n+1 normal equations.

Step 5: Solution of Normal Equations

Solve the system of normal equations using suitable techniques such as substitution, matrix methods, or Gaussian elimination to obtain the values of the coefficients.

Step 6: Construction of Regression Polynomial

Substitute the computed coefficients into the polynomial equation to obtain the final regression model. This polynomial can be used to estimate or predict values of 𝑦 for given values of x within or slightly beyond the range of the data.

Advantages

1.Can model non-linear relationships effectively

2.Flexible and adaptable to various data trends

3.Useful for interpolation and prediction

Limitations

1.High-degree polynomials may lead to overfitting

2.Poor extrapolation outside the data range

3.Sensitive to outliers


Applications

1.Engineering Analysis – Modeling stress–strain relationships and system behavior

2.Physics – Curve fitting of experimental data

3.Economics – Trend analysis and forecasting

4.Data Science – Pattern recognition and predictive modeling

5.Scientific Research – Empirical data approximation

#### Polynomial Equation Code
```cpp
#include <iostream>
#include <vector>
#include <cmath>
#include<fstream>
using namespace std;

vector<double> gaussElimination(vector<vector<double>> A, vector<double> B, ifstream & in, ofstream &out)
{
    int n = B.size();

    for(int i=0;i<n;i++)
    {
        if(fabs(A[i][i]) < 1e-9)
        {
            for(int k=i+1;k<n;k++)
            {
                if(fabs(A[k][i]) > 1e-9)
                {
                    swap(A[i], A[k]);
                    swap(B[i], B[k]);
                    break;
                }
            }
        }

        for(int k=i+1;k<n;k++)
        {
            double factor = A[k][i] / A[i][i];
            for(int j=i;j<n;j++)
                A[k][j] -= factor * A[i][j];

            B[k] -= factor * B[i];
        }
    }

    vector<double> X(n);
    for(int i=n-1;i>=0;i--)
    {
        X[i] = B[i];
        for(int j=i+1;j<n;j++)
            X[i] -= A[i][j] * X[j];

        X[i] /= A[i][i];
    }

    return X;
}

int main()
{
    string inputFile="Polynomial_input.txt";
    string outputFile="Polynomial_output.txt";

    ifstream in(inputFile);
    if(!in)
    {
        cout<<"Input file error!"<<endl;
        return 1;
    }

    ofstream out(outputFile);
    if(!out)
    {
         cout<<"Output file error!"<<endl;
        return 1;
    }
    int n, degree;
   // "Enter number of data points: ";
    in >> n;

    vector<double> x(n), y(n);

 //"Enter x values:\n";
    for(int i=0;i<n;i++)
        in >> x[i];

// "Enter y values:\n";
    for(int i=0;i<n;i++)
        in >> y[i];

    // "Enter degree of polynomial: ";
    in >> degree;

    int m = degree;
    vector<vector<double>> A(m+1, vector<double>(m+1));
    vector<double> B(m+1);

    for(int i=0;i<=m;i++)
    {
        for(int j=0;j<=m;j++)
        {
            A[i][j] = 0;
            for(int k=0;k<n;k++)
                A[i][j] += pow(x[k], i+j);
        }

        B[i] = 0;
        for(int k=0;k<n;k++)
            B[i] += pow(x[k], i) * y[k];
    }

    vector<double> coeff = gaussElimination(A, B,in,out);

 out<<"Fitted Polynomial:\n";
    out << "y = ";
    for(int i=0;i<=m;i++)
    {
        out << coeff[i];
        if(i > 0) out << "x^" << i;
        if(i != m) out << " + ";
    }
    out << endl;
   in.close();
out.close();
    return 0;
}

```

#### Polynomial Equation Input
```
5
0 1 2 3 4
1 2 5 10 17
2
```

#### Polynomial Equation Output
```
Fitted Polynomial:
y = 1 + 0x^1 + 1x^2
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
```cpp
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
Newton’s Forward Interpolation method is a numerical technique used to estimate the value of a function when the independent variable values are equally spaced and the required value lies near the beginning of the data table. This method is based on forward finite differences, which represent the successive changes in function values as the independent variable increases. By constructing a forward difference table and using Newton’s forward interpolation polynomial, the function can be approximated accurately near the starting point of the data.

Mathematical Formulation

Let the function values be given at equally spaced points: x0,x1,x2,…,xnwith corresponding values: y0,y1,y2,…,yn
	​
The step size is:h=x1−x0
	​
Define: u=(x−x0)/h
	​
The Newton’s Forward Interpolation Formula is:

y=y0+uΔy0+ 1/2!u(u−1)Δ2y0+ 1/3!u(u−1)(u−2)Δ3y0+⋯ where Δ denotes forward differences.

Procedure (Newton’s Forward Interpolation):

Step 1: Arrange the Given Data

Write the given values of 𝑥 in ascending order: x0,x1,x2,…,xn
	​
Write the corresponding values of y:y0,y1,y2,…,yn
	​

Step 2: Check Equal Spacing

Verify that: x1−x0=x2−x1= ⋯ =h . If spacing is not equal, do not use this method.

Step 3: Construct the Forward Difference Table

Compute first forward differences: Δy0=y1−y0,Δy1=y2−y1,…

Compute second forward differences: Δ2y0=Δy1−Δy0
	​
Continue until the required order of differences is obtained. Tabulate all differences neatly.

Step 4: Compute the Parameter u and Identify the value of x at which interpolation is required.​

Step 5: Write the Newton Forward Formula

Step 6: Substitute and Calculate

Substitute values of u and the differences.Use as many terms as required for accuracy.Perform calculations step by step to obtain the interpolated value. Evaluate the polynomial to obtain the required interpolated value.
#### Newton's Forward Interpolation Code
```cpp
#include <iostream>
#include<bits/stdc++.h>
#include<fstream>
using namespace std;

int main()

{

 string inputFile="Forward_equal_input.txt";
    string outputFile="Forward_equal_output.txt";

    ifstream in(inputFile);
    if(!in)
    {
        cout<<"Input file error!"<<endl;
        return 1;
    }

    ofstream out(outputFile);
    if(!out)
    {
         cout<<"Output file error!"<<endl;
        return 1;
    }
    int n;
    //"How many numbers:";
    in>>n;
    vector<double>x(n);
    vector<double>y(n);

  //"value for x: ";
    for(int i=0;i<n;i++)
        in>>x[i];

    double diff=x[1]-x[0];

    for(int i=2;i<n;i++)
    {
        if((x[i]-x[i-1])!= diff)
           {
               out<<"Difference are not equal\n";
               return 0;
           }
    }

   //"value for y: ";
    for(int i=0;i<n;i++)
        in>>y[i];


    vector<vector<double>>F(n,vector<double>(n,0));

    for(int i=0;i<n;i++)
        F[i][0]=y[i];

    //difference table calculation

    for(int j=1;j<n;j++)
    {
        for(int i=0;i<n-j;i++)
            F[i][j]=F[i+1][j-1]-F[i][j-1];
    }

    out<<"Difference table:\n";

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
            out<<F[i][j]<<"   ";
        out<<endl;
    }


    double result=0;
    double X;
   //"Give a value of x: ";
    in>>X;
    double h= x[1]-x[0];
    double u=(X-x[0])/h;

    for(int i=0;i<n;i++)
    {
        double key=1;
        for(int j=0;j<i;j++)
        {
            key*=(u-j)/(j+1);
        }

        result+=key*F[0][i];
    }

    out<<"\nResult:  "<<result<<endl;
     in.close();
out.close();

    return 0;
}

```
#### Newton's Forward Interpolation Input
```
4
1 3 5 7
4 7 9 18
2

```
#### Newton's Forward Interpolation Output
```
Difference table:
4   3   -1   8   
7   2   7   0   
9   9   0   0   
18   0   0   0   

Result:  6.125

```

### Newton's Backward Interpolation
#### Newton's Backward Interpolation Theory

Newton’s Backward Interpolation method is a numerical technique used when the data points are equally spaced and the required value lies near the end of the data table. Similar to forward interpolation, this method is based on finite differences; however, it uses backward differences, which consider changes in function values in the reverse direction.

This method is particularly useful when interpolation is required close to the last given data point, where forward differences may lead to reduced accuracy. By constructing a backward difference table and applying Newton’s backward interpolation polynomial, accurate estimates of function values can be obtained near the end of the interval. Newton’s backward interpolation is commonly used in engineering and scientific computations involving tabulated data.

Step-by-Step Procedure
Step 1: Arrange the Data

Write the x-values in ascending order.

Write corresponding y-values.

Step 2: Verify Equal Spacing

Ensure: x1−x0=x2−x1=h

Step 3: Construct the Backward Difference Table

Compute first backward differences:
∇𝑦𝑛=𝑦𝑛−𝑦𝑛−1     ∇yn=yn−yn−1
	​
Compute second backward differences:
∇𝑦𝑛=∇𝑦𝑛−∇𝑦𝑛−1       ∇2yn=∇yn−∇yn−1  Continue for higher-order differences.

Step 4: Compute the Parameter : 
v Identify the value of x where interpolation is required.

Calculate:vu=(x−xn)/h	​

Step 5: Write the Newton Backward Formula
y=yn+ v∇yn+ 1/2!v(v+1)∇2yn+ 1/3!v(v+1)(v+2)∇3yn+⋯

Step 6: Substitute and Calculate:

Substitute the values of u and backward differences.

Evaluate term by term.Obtain the required interpolated value.
#### Newton's Backward Interpolation Code
```cpp
#include <iostream>
#include<bits/stdc++.h>
#include<fstream>

using namespace std;

int main()
{
     string inputFile="Backward_input.txt";
    string outputFile="Backward_output.txt";

    ifstream in(inputFile);
    if(!in)
    {
        cout<<"Input file error!"<<endl;
        return 1;
    }

    ofstream out(outputFile);
    if(!out)
    {
         cout<<"Output file error!"<<endl;
        return 1;
    }

    int n;
    //"How many numbers:";
    in>>n;
    vector<double>x(n);
    vector<double>y(n);

 //"value for x: ";
    for(int i=0;i<n;i++)
        in>>x[i];

    double diff=x[1]-x[0];

    for(int i=2;i<n;i++)
    {
        if((x[i]-x[i-1])!= diff)
           {
               out<<"Difference are not equal\n";
               return 0;
           }
    }

    //"value for y: ";
    for(int i=0;i<n;i++)
        in>>y[i];


    vector<vector<double>>F(n,vector<double>(n,0));

    for(int i=0;i<n;i++)
        F[i][0]=y[i];

    //difference table calculation

    for(int j=1;j<n;j++)
    {
        for(int i=n-1;i>=j;i--)
            F[i][j]=F[i][j-1]-F[i-1][j-1];
    }

    out<<"Difference table:\n";

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
            out<<F[i][j]<<"   ";
        out<<endl;
    }


    double result=0;
    double X;
    //"Give a value of x: ";
    in>>X;
    double h= x[1]-x[0];
    double v=(X-x[n-1])/h;

    for(int i=0;i<n;i++)
    {
        double key=1;
        for(int j=0;j<i;j++)
        {
            key*=(v+j)/(j+1);
        }

        result+=key*F[n-1][i];
    }

    out<<"\nResult:  "<<result<<endl;
    in.close();
out.close();

    return 0;
}

```
#### Newton's Backward Interpolation Input
```
4
1 3 5 7
4 7 9 20
6

```
#### Newton's Backward Interpolation Output
```
Difference table:
4   0   0   0   
7   3   0   0   
9   2   -1   0   
20   11   9   10   

Result:  12.75

```

### Newton's Divided Difference Interpolation
#### Newton's Divided Difference Interpolation Theory

Newton’s Divided Difference Interpolation method is used when the given data points are not equally spaced. In many real-world applications, measurements are taken at irregular intervals due to experimental constraints, making finite difference methods unsuitable. Divided difference interpolation overcomes this limitation by using differences that depend directly on the spacing between data points.

This method constructs an interpolation polynomial using divided differences, which are calculated based on the ratio of differences in function values to differences in independent variable values. A major advantage of this method is its flexibility—it can be applied to both equal and unequal interval data. Additionally, new data points can be added without recalculating the entire interpolation polynomial, making it efficient and adaptable for practical data analysis.

Step-by-Step Procedure

Step 1: Arrange the Data

Write the given data points in ascending order of x:(x0,y0),(x1,y1),…,(xn,yn)

Step 2: Construct the Divided Difference Table

Compute first divided differences:

𝑓[𝑥0,𝑥]=𝑦1−𝑦0/𝑥1−𝑥0
	
Compute second divided differences:

𝑓[𝑥0,𝑥1,𝑥2]=𝑓[𝑥1,𝑥2]−𝑓[𝑥0,𝑥1]/𝑥2−𝑥0
	​Continue until all divided differences are obtained.

Record values in tabular form.

Step 3: Write the Newton Divided Difference Formula

  y=y0+(x−x0)f[x0,x1]+(x−x0)(x−x1)f[x0,x1,x2]+⋯
  
Step 4: Substitute the Required Value of 𝑥

Replace x with the given interpolation point.

Substitute the divided difference values.

Step 5: Perform Calculations

Multiply terms carefully in sequence. Stop when sufficient accuracy is achieved. Obtain the interpolated value.

#### Newton's Divided Difference Interpolation Code
```cpp
#include <iostream>
#include<bits/stdc++.h>
#include<fstream>

using namespace std;

int main()
{
    string inputFile="Forward_unequal_input.txt";
    string outputFile="Forward_unequal_output.txt";

    ifstream in(inputFile);
    if(!in)
    {
        cout<<"Input file error!"<<endl;
        return 1;
    }

    ofstream out(outputFile);
    if(!out)
    {
         cout<<"Output file error!"<<endl;
        return 1;
    }
    int n;
    //"How many numbers:";
    in>>n;
    vector<double>x(n);
    vector<double>y(n);

  //"value for x: ";
    for(int i=0;i<n;i++)
        in>>x[i];

   //"value for y: ";
    for(int i=0;i<n;i++)
        in>>y[i];

    vector<vector<double>>F(n,vector<double>(n,0));

    for(int i=0;i<n;i++)
        F[i][0]=y[i];

    //difference table calculation

    for(int j=1;j<n;j++)
    {
        for(int i=0;i<n-j;i++)
            F[i][j]=(F[i+1][j-1]-F[i][j-1])/(x[i+j]-x[i]);
    }

    out<<"Difference table:\n";

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
            out<<F[i][j]<<"   ";
        out<<endl;
    }


    double result=0;
    double X;
 //"Give a value of x: ";
    in>>X;

    for(int i=0;i<n;i++)
    {
        double key=1;
        for(int j=0;j<i;j++)
        {
            key*=(X-x[j]);
        }

        result+=key*F[0][i];
    }

    out<<"\nResult:  "<<result<<endl;
     in.close();
out.close();

    return 0;
}

```
#### Newton's Divided Difference Interpolation Input
```
4
3 5 9 14
5 8 12 20
4

```
#### Newton's Divided Difference Interpolation Output
```
Difference table:
5   1.5   -0.0833333   0.0136364   
8   1   0.0666667   0   
12   1.6   0   0   
20   0   0   0   

Result:  6.65152

```

### Solution of Numerical Integration
### Simpson's One-Third Rule
#### Simpson's One-Third Rule Theory
Simpson’s 1/3 rule is a numerical integration technique that approximates the curve by a series of parabolic arcs. The rule requires the interval 
[a,b] to be divided into an even number of equal sub-intervals. Function values are evaluated at equally spaced points, and a weighted sum of these values is used to approximate the area under the curve.

Mathematical Formulation

Let
𝑦=𝑓(𝑥)  be the given function to be integrated over the closed interval [a,b].
Let the total number of sub-intervals be 𝑛, where n is even. Then the number of data points is n+1.

The step size is: ℎ=(b−a)/n
	​Let:
𝑦0=𝑓(𝑥0),  𝑦1=𝑓(𝑥1),…,𝑦𝑛=𝑓(𝑥𝑛)

The Simpson’s 1/3 rule formula is:
∫abydx ≈ 3h [y0+yn+4(y1+y3+⋯+yn−1)+2(y2+y4+⋯+yn−2)]

Procedure (Simpson’s 1/3 Rule)

1.Select n+1 equally spaced values of x in the interval [a,b], where n is even.
2.Compute the corresponding values of 𝑦 using the given function y=f(x).
3.Calculate the step size h=(b−a)/n.
4.Substitute the values of y0,y1,…,yn into Simpson’s 1/3 formula.
5.Perform the calculations to obtain the approximate value of the integral.

#### Simpson's One-Third Rule Code
```cpp
#include <iostream>
#include<bits/stdc++.h>
#include<fstream>

using namespace std;

double function_call(double x)
{
    return 1/(1+(x*x)); // for function 1/(1+ x^2)
}

double odd(vector<double>&y, int n)
{
    double sum=0;
    for(int i=1;i<n;i+=2)
        sum+=y[i];
    return sum;
}

double even(vector<double>&y ,int n)
{
    double sum=0;
    for(int i=2;i<n-1;i+=2)
        sum+=y[i];
    return sum;
}
int main()
{
     string inputFile="Simpson(1-3)_input.txt";
    string outputFile="Simpson(1-3)_output.txt";

    ifstream in(inputFile);
    if(!in)
    {
        cout<<"Input file error!"<<endl;
        return 1;
    }

    ofstream out(outputFile);
    if(!out)
    {
         cout<<"Output file error!"<<endl;
        return 1;
    }
 //"Print how many interval: ";
    int n;
    in>>n;
    if(n>=2 && n%2==0)
    {
       vector<double>x(n+1);
    vector<double>y(n+1);
    for(int i=0;i<=n;i++)
        in>>x[i];

        double check=x[1]-x[0];

      //checking if the intervals are equal or not?
        for(int i=2;i<=n;i++)
        {
            if((x[i]-x[i-1]!=check))
               {
                   out<<"Intervals are not same!\n";
                   return 0;
               }
        }

    for(int i=0;i<=n;i++)
        y[i]=function_call(x[i]);

        double a,b;
       //"Print upper-bound and lower-bound: ";
        in>>b>>a;

    double h= (b-a)/n;
    double c=y[0]+y[n];
    double ans= (h*(c + (4*odd(y,n)) + (2*even(y,n))))/3;

    out<<fixed<<setprecision(5);
    out<<"Result is: "<<ans<<endl;
    }

    else
        out<<"Interval is not even!\n";
 in.close();
out.close();
    return 0;
}

```
#### Simpson's One-Third Rule Input
```
4
3 5 7 9 11
2 5

```
#### Simpson's One-Third Rule Output
```
Result is: -0.08771

```
### Simpson's Three-Eighth Rule
#### Simpson's Three-Eighth Rule Theory
 Simpson’s 3/8 rule is another closed Newton–Cotes integration method that approximates the function using cubic polynomials. The interval is divided into groups of three equal sub-intervals, and a weighted average of function values is used to compute the area.
This method requires the total number of sub-intervals to be a multiple of 3. It is often used in combination with Simpson’s 1/3 rule when the number of intervals is not even.

Mathematical Formulation:

Let the interval [a,b] be divided into 𝑛 sub-intervals, where n is a multiple of 3.

The step size is: h= (b−a)/n

The Simpson’s 3/8 rule formula is:
∫ab ydx≈ 3h/8[y0+yn+3(y1+y2+y4+y5+⋯+yn−1)+2(y3+y6+⋯+yn−3)]

Procedure (Simpson’s 3/8 Rule)

1.Divide the interval [a,b] into n equal sub-intervals, where n is a multiple of 3.

2.Calculate the step size h=(b−a)/n.

3.Evaluate the function values y=f(x) at each of the n+1 points.

4.Substitute the values into Simpson’s 3/8 formula.

5.Compute the final approximate value of the integral.

#### Simpson's Three-Eighth Rule Code
```cpp
#include <iostream>
#include<bits/stdc++.h>
#include<fstream>

using namespace std;

double function_call(double x)
{
    return sqrt(x); // for function 1/(1+ x^2)
}

double sum1(vector<double>&y, int n)
{
    double sum=0;
    for(int i=1;i<n;i++)
    {
        if((i%3)!=0)
             sum+=y[i];
    }

    return sum;
}

double sum2(vector<double>&y ,int n)
{
    double sum=0;
    for(int i=3;i<n;i+=3)
        sum+=y[i];
    return sum;
}
int main()
{
     string inputFile="Simpson(3-8)_input.txt";
    string outputFile="Simpson(3-8)_output.txt";

    ifstream in(inputFile);
    if(!in)
    {
        cout<<"Input file error!"<<endl;
        return 1;
    }

    ofstream out(outputFile);
    if(!out)
    {
         cout<<"Output file error!"<<endl;
        return 1;
    }
   //"Print how many interval: ";
    int n;
    in>>n;
    if(n>=3)
    {
       vector<double>x(n+1);
    vector<double>y(n+1);
    for(int i=0;i<=n;i++)
        in>>x[i];

        double check=x[1]-x[0];

      //checking if the intervals are equal or not?
        for(int i=2;i<=n;i++)
        {
            if((x[i]-x[i-1]!=check))
               {
                   out<<"Intervals are not same!\n";
                   return 0;
               }
        }

    for(int i=0;i<=n;i++)
        y[i]=function_call(x[i]);

        double a,b;
        //"Print upper-bound and lower-bound: ";
        in>>b>>a;

    double h= (b-a)/n;
    double c=y[0]+y[n];
    double ans= (3*h*(c + (3*sum1(y,n)) + (2*sum2(y,n))))/8;

    for(int i=0;i<=n;i++)
        out<<y[i]<<"  ";

    out<<endl;

    out<<fixed<<setprecision(5);
    out<<"Result is: "<<ans<<endl;
    }

    else
        out<<"Interval is less than 3!\n";
  in.close();
out.close();
    return 0;
}

```
#### Simpson's Three-Eighth Rule Input
```
3
0 1 2 3
3 0

```
#### Simpson's Three-Eighth Rule Output
```
0  1  1.41421  1.73205  
Result is: 3.36551

```



---





