# include <cmath>
# include <complex>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <fstream>
# include <iomanip>
# include <iostream>

using namespace std;

# include "r8lib.hpp"



void gamma_values ( int &n_data, double &x, double &fx )































































{
# define N_MAX 25

  static double fx_vec[N_MAX] = {
     -0.3544907701811032E+01,
     -0.1005871979644108E+03,
      0.9943258511915060E+02,
      0.9513507698668732E+01,
      0.4590843711998803E+01,
      0.2218159543757688E+01,
      0.1772453850905516E+01,
      0.1489192248812817E+01,
      0.1164229713725303E+01,
      0.1000000000000000E+01,
      0.9513507698668732E+00,
      0.9181687423997606E+00,
      0.8974706963062772E+00,
      0.8872638175030753E+00,
      0.8862269254527580E+00,
      0.8935153492876903E+00,
      0.9086387328532904E+00,
      0.9313837709802427E+00,
      0.9617658319073874E+00,
      0.1000000000000000E+01,
      0.2000000000000000E+01,
      0.6000000000000000E+01,
      0.3628800000000000E+06,
      0.1216451004088320E+18,
      0.8841761993739702E+31 };

  static double x_vec[N_MAX] = {
     -0.50E+00,
     -0.01E+00,
      0.01E+00,
      0.10E+00,
      0.20E+00,
      0.40E+00,
      0.50E+00,
      0.60E+00,
      0.80E+00,
      1.00E+00,
      1.10E+00,
      1.20E+00,
      1.30E+00,
      1.40E+00,
      1.50E+00,
      1.60E+00,
      1.70E+00,
      1.80E+00,
      1.90E+00,
      2.00E+00,
      3.00E+00,
      4.00E+00,
     10.00E+00,
     20.00E+00,
     30.00E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}


void gamma_log_values ( int &n_data, double &x, double &fx )



















































{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.1524063822430784E+01,
      0.7966778177017837E+00,
      0.3982338580692348E+00,
      0.1520596783998375E+00,
      0.0000000000000000E+00,
     -0.4987244125983972E-01,
     -0.8537409000331584E-01,
     -0.1081748095078604E+00,
     -0.1196129141723712E+00,
     -0.1207822376352452E+00,
     -0.1125917656967557E+00,
     -0.9580769740706586E-01,
     -0.7108387291437216E-01,
     -0.3898427592308333E-01,
     0.00000000000000000E+00,
     0.69314718055994530E+00,
     0.17917594692280550E+01,
     0.12801827480081469E+02,
     0.39339884187199494E+02,
     0.71257038967168009E+02 };

  static double x_vec[N_MAX] = {
      0.20E+00,
      0.40E+00,
      0.60E+00,
      0.80E+00,
      1.00E+00,
      1.10E+00,
      1.20E+00,
      1.30E+00,
      1.40E+00,
      1.50E+00,
      1.60E+00,
      1.70E+00,
      1.80E+00,
      1.90E+00,
      2.00E+00,
      3.00E+00,
      4.00E+00,
     10.00E+00,
     20.00E+00,
     30.00E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}


int i4_log_10 ( int i )

















































{
  int i_abs;
  int ten_pow;
  int value;

  if ( i == 0 )
  {
    value = 0;
  }
  else
  {
    value = 0;
    ten_pow = 10;

    i_abs = abs ( i );

    while ( ten_pow <= i_abs )
    {
      value = value + 1;
      ten_pow = ten_pow * 10;
    }

  }

  return value;
}


int i4_max ( int i1, int i2 )

























{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}


int i4_min ( int i1, int i2 )

























{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}


int i4_modp ( int i, int j )




















































{
  int value;

  if ( j == 0 )
  {
    cerr << "\n";
    cerr << "I4_MODP - Fatal error!\n";
    cerr << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}


int i4_power ( int i, int j )

























{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}


int i4_sign ( int i )
























{
  int value;

  if ( i < 0 )
  {
    value = -1;
  }
  else
  {
    value = 1;
  }
  return value;
}


int i4_uniform_ab ( int a, int b, int &seed )






























































{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }



  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;



  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );



  value = round ( r );



  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}


int i4_wrap ( int ival, int ilo, int ihi )



















































{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i4_min ( ilo, ihi );
  jhi = i4_max ( ilo, ihi );

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

  return value;
}


double i4int_to_r8int ( int imin, int imax, int i, double rmin, double rmax )



































{
  double r;

  if ( imax == imin )
  {
    r = 0.5 * ( rmin + rmax );
  }
  else
  {
    r = ( ( double ) ( imax - i        ) * rmin
        + ( double ) (        i - imin ) * rmax )
        / ( double ) ( imax     - imin );
  }

  return r;
}


void i4mat_print ( int m, int n, int a[], string title )

































{
  i4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}


void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, string title )






































{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }



  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    cout << "\n";





    cout << "  Col:";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << "  " << setw(6) << j - 1;
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";



    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( ihi < m )
    {
      i2hi = ihi;
    }
    else
    {
      i2hi = m;
    }


    for ( i = i2lo; i <= i2hi; i++ )
    {



      cout << setw(5) << i - 1 << ":";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << "  " << setw(6) << a[i-1+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}


void i4vec_copy ( int n, int a1[], int a2[] )































{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}


int *i4vec_indicator0_new ( int n )





























{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = i;
  }
  return a;
}


int *i4vec_indicator1_new ( int n )





























{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = i + 1;
  }
  return a;
}


void i4vec_permute ( int n, int p[], int a[] )




















































{
  int a_temp;
  int i;
  int iget;
  int iput;
  int istart;

  if ( !perm0_check ( n, p ) )
  {
    cerr << "\n";
    cerr << "I4VEC_PERMUTE - Fatal error!\n";
    cerr << "  PERM0_CHECK rejects permutation.\n";
    exit ( 1 );
  }





  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1;
  }



  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      a_temp = a[istart-1];
      iget = istart;



      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cerr << "\n";
          cerr << "I4VEC_PERMUTE - Fatal error!\n";
          cerr << "  Entry IPUT = " << iput << " of the permutation has\n";
          cerr << "  an illegal value IGET = " << iget << ".\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[iput-1] = a_temp;
          break;
        }
        a[iput-1] = a[iget-1];
      }
    }
  }



  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }



  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1;
  }

  return;
}


void i4vec_print ( int n, int a[], string title )































{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8) << i
         << ": " << setw(8) << a[i]  << "\n";
  }
  return;
}


void i4vec_transpose_print ( int n, int a[], string title )








































{
  int i;
  int ihi;
  int ilo;
  int title_len;

  title_len = title.length ( );

  for ( ilo = 1; ilo <= n; ilo = ilo + 5 )
  {
    ihi = ilo + 5 - 1;
    if ( n < ihi )
    {
      ihi = n;
    }

    if ( ilo == 1 )
    {
      cout << title;
    }
    else
    {
      for ( i = 1; i <= title_len; i++ )
      {
        cout << " ";
      }
    }
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << setw(12) << a[i-1];
    }
    cout << "\n";
  }

  return;
}


void i4vec_zeros ( int n, int a[] )





























{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return;
}


int *i4vec_zeros_new ( int n )





























{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return a;
}


double *legendre_zeros ( int order )




































{
  double d1;
  double d2pn;
  double d3pn;
  double d4pn;
  double dp;
  double dpn;
  double e1;
  double fx;
  double h;
  int i;
  int iback;
  int k;
  int m;
  int mp1mi;
  int ncopy;
  int nmove;
  double p;
  double pk;
  double pkm1;
  double pkp1;
  const double r8_pi = 3.141592653589793;
  double t;
  double u;
  double v;
  double x0;
  double *xtab;
  double xtemp;

  xtab = new double[order];

  e1 = ( double ) ( order * ( order + 1 ) );

  m = ( order + 1 ) / 2;

  for ( i = 1; i <= m; i++ )
  {
    mp1mi = m + 1 - i;

    t = ( double ) ( 4 * i - 1 ) * r8_pi / ( double ) ( 4 * order + 2 );

    x0 = cos ( t ) * ( 1.0 - ( 1.0 - 1.0 / ( double ) ( order ) ) 
      / ( double ) ( 8 * order * order ) );

    pkm1 = 1.0;
    pk = x0;

    for ( k = 2; k <= order; k++ )
    {
      pkp1 = 2.0 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) / ( double ) ( k );
      pkm1 = pk;
      pk = pkp1;
    }

    d1 = ( double ) ( order ) * ( pkm1 - x0 * pk );

    dpn = d1 / ( 1.0 - x0 * x0 );

    d2pn = ( 2.0 * x0 * dpn - e1 * pk ) / ( 1.0 - x0 * x0 );

    d3pn = ( 4.0 * x0 * d2pn + ( 2.0 - e1 ) * dpn ) / ( 1.0 - x0 * x0 );

    d4pn = ( 6.0 * x0 * d3pn + ( 6.0 - e1 ) * d2pn ) / ( 1.0 - x0 * x0 );

    u = pk / dpn;
    v = d2pn / dpn;



    h = -u * ( 1.0 + 0.5 * u * ( v + u * ( v * v - d3pn / ( 3.0 * dpn ) ) ) );



    p = pk + h * ( dpn + 0.5 * h * ( d2pn + h / 3.0 
      * ( d3pn + 0.25 * h * d4pn ) ) );

    dp = dpn + h * ( d2pn + 0.5 * h * ( d3pn + h * d4pn / 3.0 ) );

    h = h - p / dp;

    xtemp = x0 + h;

    xtab[mp1mi-1] = xtemp;

    fx = d1 - h * e1 * ( pk + 0.5 * h * ( dpn + h / 3.0 
      * ( d2pn + 0.25 * h * ( d3pn + 0.2 * h * d4pn ) ) ) );
  }

  if ( ( order % 2 ) == 1 )
  {
    xtab[0] = 0.0;
  }



  nmove = ( order + 1 ) / 2;
  ncopy = order - nmove;

  for ( i = 1; i <= nmove; i++ )
  {
    iback = order + 1 - i;
    xtab[iback-1] = xtab[iback-ncopy-1];
  }



  for ( i = 1; i <= order - nmove; i++ )
  {
    xtab[i-1] = - xtab[order-i];
  }

  return xtab;
}


bool perm0_check ( int n, int p[] )


































{
  bool check;
  int location;
  int value;

  check = true;

  for ( value = 0; value < n; value++ )
  {
    check = false;

    for ( location = 0; location < n; location++ )
    {
      if ( p[location] == value )
      {
        check = true;
        break;
      }
    }

    if ( ! check )
    {
      cout << "\n";
      cout << "PERM0_CHECK - Fatal error!\n";
      cout << "  Permutation is missing value " << value << "\n";
      break;
    }

  }

  return check;
}


int *perm0_uniform_new ( int n, int &seed )



































{
  int i;
  int j;
  int k;
  int *p;

  p = new int[n];

  for ( i = 0; i < n; i++ )
  {
    p[i] = i;
  }

  for ( i = 0; i < n; i++ )
  {
    j = i4_uniform_ab ( i, n - 1, seed );
    k    = p[i];
    p[i] = p[j];
    p[j] = k;
  }

  return p;
}


bool perm1_check ( int n, int p[] )


































{
  bool check;
  int location;
  int value;

  check = true;

  for ( value = 1; value <= n; value++ )
  {
    check = false;

    for ( location = 0; location < n; location++ )
    {
      if ( p[location] == value )
      {
        check = true;
        break;
      }
    }

    if ( ! check )
    {
      cout << "\n";
      cout << "PERM1_CHECK - Fatal error!\n";
      cout << "  Permutation is missing value " << value << "\n";
      break;
    }

  }

  return check;
}


int *perm1_uniform_new ( int n, int &seed )



































{
  int i;
  int j;
  int k;
  int *p;

  p = new int[n];

  for ( i = 0; i < n; i++ )
  {
    p[i] = i + 1;
  }

  for ( i = 0; i < n; i++ )
  {
    j = i4_uniform_ab ( i, n - 1, seed );
    k    = p[i];
    p[i] = p[j];
    p[j] = k;
  }

  return p;
}


double r8_abs ( double x )





























{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}


double r8_acos ( double c )































{
  const double r8_pi = 3.141592653589793;
  double value;

  if ( c <= -1.0 )
  {
    value = r8_pi;
  }
  else if ( 1.0 <= c )
  {
    value = 0.0;
  }
  else
  {
    value = acos ( c );
  }
  return value;
}


double r8_acosh ( double x )



























































{
  double value;

  if ( x < 1.0 )
  {
    cerr << "\n";
    cerr << "R8_ACOSH - Fatal error!\n";
    cerr << "  Argument X must satisfy 1 <= X.\n";
    cerr << "  The input X = " << x << "\n";
    exit ( 1 );
  }

  value = 2.0 * log ( 
    sqrt ( 0.5 * ( x + 1.0 ) ) + sqrt ( 0.5 * ( x - 1.0 ) ) );

  return value;
}


double r8_add ( double x, double y )

























{
  double value;

  value = x + y;

  return value;
}


double r8_agm ( double a, double b )





















































{
  double a1;
  double a2;
  double b1;
  double b2;
  int it;
  int it_max = 1000;
  double tol;
  double value;

  if ( a < 0.0 )
  {
    cerr << "\n";
    cerr << "R8_AGM - Fatal error!\n";
    cerr << "  A < 0.\n";
    exit ( 1 );
  }

  if ( b < 0.0 )
  {
    cerr << "\n";
    cerr << "R8_AGM - Fatal error!\n";
    cerr << "  B < 0.\n";
    exit ( 1 );
  }

  if ( a == 0.0 || b == 0.0 )
  {
    value = 0.0;
    return value;
  }

  if ( a == b )
  {
    value = a;
    return value;
  }

  a1 = a;
  b1 = b;

  it = 0;
  tol = 100.0 * r8_epsilon ( );

  for ( ; ; )
  {
    it = it + 1;

    a2 = ( a1 + b1 ) / 2.0;
    b2 = sqrt ( a1 * b1 );

    if ( fabs ( a2 - b2 ) <= tol * ( a2 + b2 ) )
    {
      break;
    }

    if ( it_max < it )
    {
      break;
    }

    a1 = a2;
    b1 = b2;
  }
  value = a2;

  return value;
}


double r8_aint ( double x )

























{
  double value;

  if ( x < 0.0 )
  {
    value = - ( double ) ( ( int ) ( fabs ( x ) ) );
  }
  else
  {
    value =   ( double ) ( ( int ) ( fabs ( x ) ) );
  }

  return value;
}


double r8_asin ( double s )































{
  double angle;
  const double r8_pi = 3.141592653589793;

  if ( s <= -1.0 )
  {
    angle = - r8_pi / 2.0;
  }
  else if ( 1.0 <= s )
  {
    angle = r8_pi / 2.0;
  }
  else
  {
    angle = asin ( s );
  }
  return angle;
}


double r8_asinh ( double x )




































{
  double value;

  value = log ( x + sqrt ( x * x + 1.0 ) );

  return value;
}


double r8_atan ( double y, double x )











































{
  double abs_x;
  double abs_y;
  const double r8_pi = 3.141592653589793;
  double theta = 0.0;
  double theta_0;



  if ( x == 0.0 )
  {
    if ( 0.0 < y )
    {
      theta = r8_pi / 2.0;
    }
    else if ( y < 0.0 )
    {
      theta = 3.0 * r8_pi / 2.0;
    }
    else if ( y == 0.0 )
    {
      theta = 0.0;
    }
  }
  else if ( y == 0.0 )
  {
    if ( 0.0 < x )
    {
      theta = 0.0;
    }
    else if ( x < 0.0 )
    {
      theta = r8_pi;
    }
  }



  else
  {
    abs_y = fabs ( y );
    abs_x = fabs ( x );

    theta_0 = atan2 ( abs_y, abs_x );

    if ( 0.0 < x && 0.0 < y )
    {
      theta = theta_0;
    }
    else if ( x < 0.0 && 0.0 < y )
    {
      theta = r8_pi - theta_0;
    }
    else if ( x < 0.0 && y < 0.0 )
    {
      theta = r8_pi + theta_0;
    }
    else if ( 0.0 < x && y < 0.0 )
    {
      theta = 2.0 * r8_pi - theta_0;
    }
  }

  return theta;
}


double r8_atanh ( double x )



































{
  const double r8_huge = 1.79769313486231571E+308;
  double value;

  if ( x <= -1.0 )
  {
    value = - r8_huge;
  }
  else if ( 1.0 <= x )
  {
    value = + r8_huge;
  }
  else
  {
    value = 0.5 * log ( ( 1.0 + x ) / ( 1.0 - x ) );
  }

  return value;
}


double r8_big ( )





























{
  double value;

  value = 1.0E+30;

  return value;
}


double r8_cas ( double x )



































{
  double value;

  value = cos ( x ) + sin ( x );

  return value;
}


double r8_ceiling ( double x )










































{
  double value;

  value = ( double ) ( ( int ) x );

  if ( value < x )
  {
    value = value + 1.0;
  }

  return value;
}


double r8_choose ( int n, int k )

















































{
  double arg;
  double fack;
  double facn;
  double facnmk;
  double value;

  if ( n < 0 )
  {
    value = 0.0;
  }
  else if ( k == 0 )
  {
    value = 1.0;
  }
  else if ( k == 1 )
  {
    value = ( double ) ( n );
  }
  else if ( 1 < k && k < n - 1 )
  {
    arg = ( double ) ( n + 1 );
    facn = r8_gamma_log ( arg );

    arg = ( double ) ( k + 1 );
    fack = r8_gamma_log ( arg );

    arg = ( double ) ( n - k + 1 );
    facnmk = r8_gamma_log ( arg );

    value = r8_nint ( exp ( facn - fack - facnmk ) );
  }
  else if ( k == n - 1 )
  {
    value = ( double ) ( n );
  }
  else if ( k == n )
  {
    value = 1.0;
  }
  else
  {
    value = 0.0;
  }

  return value;
}


double r8_chop ( int place, double x )














































{
  double fac;
  int temp;
  double value;

  temp = ( int ) ( r8_log_2 ( x ) );
  fac = pow ( 2.0, ( temp - place + 1 ) );
  value = ( double ) ( ( int ) ( x / fac ) ) * fac;

  return value;
}


double r8_cosd ( double degrees )

























{
  const double r8_pi = 3.141592653589793;
  double radians;
  double value;

  radians = r8_pi * ( degrees / 180.0 );

  value = cos ( radians );

  return value;
}


double r8_cot ( double angle )





























{
  double value;

  value = cos ( angle ) / sin ( angle );

  return value;
}


double r8_cotd ( double degrees )

























{
  const double r8_pi = 3.141592653589793;
  double radians;
  double value;

  radians = r8_pi * ( degrees / 180.0 );

  value = cos ( radians ) / sin ( radians );

  return value;
}


double r8_csc ( double theta )






























{
  double value;

  value = sin ( theta );

  if ( value == 0.0 )
  {
    cerr << " \n";
    cerr << "R8_CSC - Fatal error!\n";
    cerr << "  Cosecant undefined for THETA = " << theta << "\n";
    exit ( 1 );
  }

  value = 1.0 / value;

  return value;
}


double r8_cscd ( double degrees )

























{
  const double r8_pi = 3.141592653589793;
  double radians;
  double value;

  radians = r8_pi * ( degrees / 180.0 );

  value = 1.0 / sin ( radians );

  return value;
}


double r8_cube_root ( double x )






























{
  double value;

  if ( 0.0 < x )
  {
    value = pow ( ( double ) x, ( 1.0 / 3.0 ) );
  }
  else if ( x == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = - pow ( ( double ) ( fabs ( x ) ), ( 1.0 / 3.0 ) );
  }

  return value;
}


double r8_degrees ( double radians )

























{
  const double r8_pi = 3.1415926535897932384626434;
  double value;

  value = radians * 180.0 / r8_pi;

  return value;
}


double r8_diff ( double x, double y, int n )







































{
  double cx;
  double cy;
  double pow2;
  double size;
  double value;

  if ( x == y )
  {
    value = 0.0;
    return value;
  }

  pow2 = pow ( 2.0, n );




  size = r8_max ( fabs ( x ), fabs ( y ) );




  cx = x / size;
  cy = y / size;






  cx = ( double ) ( ( int ) ( cx * pow2 + 0.5 * r8_sign ( cx ) ) ) / pow2;
  cy = ( double ) ( ( int ) ( cy * pow2 + 0.5 * r8_sign ( cy ) ) ) / pow2;



  value = cx - cy;



  value = value * size;

  return value;
}


int r8_digit ( double x, int idigit )






























{
  int digit;
  int i;
  int ival;

  if ( x == 0.0 )
  {
    digit = 0;
    return digit;
  }

  if ( idigit <= 0 )
  {
    digit = 0;
    return digit;
  }



  x = fabs ( x );

  while ( x < 1.0 )
  {
    x = x * 10.0;
  }

  while ( 10.0 <= x )
  {
    x = x / 10.0;
  }

  for ( i = 1; i <= idigit; i++ )
  {
    ival = ( int ) ( x );
    x = ( x - ( double ) ival ) * 10.0;
  }

  digit = ival;

  return digit;
}


double r8_divide_i4 ( int  i, int j )

























{
  double value;

  value = ( double ) ( i ) / ( double ) ( j );

  return value;
}


double r8_e ( )



























{
  const double r8_e_save = 2.718281828459045235360287;
  double value;

  value = r8_e_save;

  return value;
}


double r8_epsilon ( )































{
  const double value = 2.220446049250313E-016;

  return value;
}


double r8_epsilon_compute ( )































{
  double one;
  double temp;
  double test;
  static double value = 0.0;

  if ( value == 0.0 )
  {
    one = ( double ) ( 1 );

    value = one;
    temp = value / 2.0;
    test = r8_add ( one, temp );

    while ( one < test )
    {
      value = temp;
      temp = value / 2.0;
      test = r8_add ( one, temp );
    }
  }

  return value;
}


double r8_exp ( double x )





































{
  const double r8_big = 1.0E+30;
  const double r8_log_max = +69.0776;
  const double r8_log_min = -69.0776;
  double value;

  if ( x <= r8_log_min )
  {
    value = 0.0;
  }
  else if ( x < r8_log_max )
  {
    value = exp ( x );
  }
  else
  {
    value = r8_big;
  }

  return value;
}


double r8_factorial ( int n )






























{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}


double r8_factorial_stirling ( int n )




































{
  const double r8_e = 2.71828182845904523;
  const double r8_pi = 3.14159265358979323;
  double value;

  if ( n < 0 )
  {
    value = 0.0;
  }
  else if ( n == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = sqrt ( 2.0 * r8_pi * ( double ) ( n ) )
      * pow ( ( double ) ( n ) / r8_e, n )
      * exp ( 1.0 / ( double ) ( 12 * n ) );
  }

  return value;
}


void r8_factorial_values ( int &n_data, int &n, double &fn )



























































{
# define N_MAX 25

  static double fn_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.6000000000000000E+01,
     0.2400000000000000E+02,
     0.1200000000000000E+03,
     0.7200000000000000E+03,
     0.5040000000000000E+04,
     0.4032000000000000E+05,
     0.3628800000000000E+06,
     0.3628800000000000E+07,
     0.3991680000000000E+08,
     0.4790016000000000E+09,
     0.6227020800000000E+10,
     0.8717829120000000E+11,
     0.1307674368000000E+13,
     0.2092278988800000E+14,
     0.3556874280960000E+15,
     0.6402373705728000E+16,
     0.1216451004088320E+18,
     0.2432902008176640E+19,
     0.1551121004333099E+26,
     0.3041409320171338E+65,
     0.9332621544394415E+158,
     0.5713383956445855E+263 };

  static int n_vec[N_MAX] = {
       0,
       1,
       2,
       3,
       4,
       5,
       6,
       7,
       8,
       9,
      10,
      11,
      12,
      13,
      14,
      15,
      16,
      17,
      18,
      19,
      20,
      25,
      50,
     100,
     150 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    fn = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    fn = fn_vec[n_data-1];
  }

  return;
# undef N_MAX
}


double r8_factorial2 ( int n )















































{
  int n_copy;
  double value;

  value = 1.0;

  if ( n < 1 )
  {
    return value;
  }

  n_copy = n;

  while ( 1 < n_copy )
  {
    value = value * ( double ) n_copy;
    n_copy = n_copy - 2;
  }

  return value;
}


void r8_factorial2_values ( int &n_data, int &n, double &f )











































































{
# define N_MAX 16

  static double f_vec[N_MAX] = {
          1.0,
          1.0,
          2.0,
          3.0,
          8.0,
         15.0,
         48.0,
        105.0,
        384.0,
        945.0,
       3840.0,
      10395.0,
      46080.0,
     135135.0,
     645120.0,
    2027025.0 };

  static int n_vec[N_MAX] = {
     0,
     1,  2,  3,  4,  5,
     6,  7,  8,  9, 10,
    11, 12, 13, 14, 15 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    f = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    f = f_vec[n_data-1];
  }

  return;
# undef N_MAX
}


double r8_fall ( double x, int n )













































{
  int i;
  double value;

  value = 1.0;

  if ( 0 < n )
  {
    for ( i = 1; i <= n; i++ )
    {
      value = value * x;
      x = x - 1.0;
    }
  }
  else if ( n < 0 )
  {
    for ( i = -1; n <= i; i-- )
    {
      value = value * x;
      x = x + 1.0;
    }
  }

  return value;
}


void r8_fall_values ( int &n_data, double &x, int &n, double &f )



















































{
# define N_MAX 15

  static double f_vec[N_MAX] = {
    120.0000000000000,
    163.1601562500000,
    216.5625000000000,
    281.6601562500000,
    360.0000000000000,
    1.000000000000000,
    7.500000000000000,
    48.75000000000000,
    268.1250000000000,
    1206.562500000000,
    4222.968750000000,
    10557.42187500000,
    15836.13281250000,
    7918.066406250000,
    -3959.03320312500 };

  static int n_vec[N_MAX] = {
    4,
    4,
    4,
    4,
    4,
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9 };

  static double x_vec[N_MAX] = {
    5.00,
    5.25,
    5.50,
    5.75,
    6.00,
    7.50,
    7.50,
    7.50,
    7.50,
    7.50,
    7.50,
    7.50,
    7.50,
    7.50,
    7.50 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    n = 0;
    f = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    n = n_vec[n_data-1];
    f = f_vec[n_data-1];
  }

  return;
# undef N_MAX
}


double r8_floor ( double x )










































{
  double value;

  value = ( double ) ( ( int ) x );

  if ( x < value )
  {
    value = value - 1.0;
  }

  return value;
}


double r8_fraction ( int i, int j )










































{
  double value;

  value = ( double ) ( i ) / ( double ) ( j );

  return value;
}


double r8_fractional ( double x )

























{
  double value;

  value = fabs ( x ) - ( double ) ( ( int ) fabs ( x ) );

  return value;
}


double r8_gamma ( double x )























































{
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double fact;
  int i;
  int n;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  bool parity;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  const double r8_epsilon = 2.220446049250313E-016;
  const double r8_pi = 3.1415926535897932384626434;
  double res;
  const double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;

  parity = false;
  fact = 1.0;
  n = 0;
  y = x;



  if ( y <= 0.0 )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != 0.0 )
    {
      if ( y1 != ( double ) ( int ) ( y1 * 0.5 ) * 2.0 )
      {
        parity = true;
      }

      fact = - r8_pi / sin ( r8_pi * res );
      y = y + 1.0;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }



  if ( y < r8_epsilon )
  {



    if ( xminin <= y )
    {
      res = 1.0 / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < 12.0 )
  {
    y1 = y;



    if ( y < 1.0 )
    {
      z = y;
      y = y + 1.0;
    }




    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - 1.0;
    }



    xnum = 0.0;
    xden = 1.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + 1.0;



    if ( y1 < y )
    {
      res = res / y1;
    }



    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + 1.0;
      }
    }
  }
  else
  {



    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - 0.5 ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }



  if ( parity )
  {
    res = - res;
  }

  if ( fact != 1.0 )
  {
    res = fact / res;
  }

  value = res;

  return value;
}


double r8_gamma_log ( double x )























































{
  double c[7] = {
    -1.910444077728E-03, 
     8.4171387781295E-04, 
    -5.952379913043012E-04,
     7.93650793500350248E-04, 
    -2.777777777777681622553E-03, 
     8.333333333333333331554247E-02, 
     5.7083835261E-03 };
  double corr;
  const double d1 = -5.772156649015328605195174E-01;
  const double d2 = 4.227843350984671393993777E-01;
  const double d4 = 1.791759469228055000094023;
  const double frtbig = 2.25E+76;
  int i;
  double p1[8] = {
    4.945235359296727046734888, 
    2.018112620856775083915565E+02, 
    2.290838373831346393026739E+03, 
    1.131967205903380828685045E+04, 
    2.855724635671635335736389E+04, 
    3.848496228443793359990269E+04, 
    2.637748787624195437963534E+04, 
    7.225813979700288197698961E+03 };
  double p2[8] = { 
    4.974607845568932035012064, 
    5.424138599891070494101986E+02, 
    1.550693864978364947665077E+04, 
    1.847932904445632425417223E+05, 
    1.088204769468828767498470E+06, 
    3.338152967987029735917223E+06, 
    5.106661678927352456275255E+06, 
    3.074109054850539556250927E+06 };
  double p4[8] = {
    1.474502166059939948905062E+04, 
    2.426813369486704502836312E+06, 
    1.214755574045093227939592E+08, 
    2.663432449630976949898078E+09, 
    2.940378956634553899906876E+10, 
    1.702665737765398868392998E+11, 
    4.926125793377430887588120E+11, 
    5.606251856223951465078242E+11 };
  double q1[8] = { 
    6.748212550303777196073036E+01, 
    1.113332393857199323513008E+03, 
    7.738757056935398733233834E+03, 
    2.763987074403340708898585E+04, 
    5.499310206226157329794414E+04, 
    6.161122180066002127833352E+04, 
    3.635127591501940507276287E+04, 
    8.785536302431013170870835E+03 };
  double q2[8] = { 
    1.830328399370592604055942E+02, 
    7.765049321445005871323047E+03, 
    1.331903827966074194402448E+05, 
    1.136705821321969608938755E+06, 
    5.267964117437946917577538E+06, 
    1.346701454311101692290052E+07, 
    1.782736530353274213975932E+07, 
    9.533095591844353613395747E+06 };
  double q4[8] = { 
    2.690530175870899333379843E+03, 
    6.393885654300092398984238E+05, 
    4.135599930241388052042842E+07, 
    1.120872109616147941376570E+09, 
    1.488613728678813811542398E+10, 
    1.016803586272438228077304E+11, 
    3.417476345507377132798597E+11, 
    4.463158187419713286462081E+11 };
  const double r8_epsilon = 2.220446049250313E-016;
  double res;
  const double sqrtpi = 0.9189385332046727417803297;
  const double xbig = 2.55E+305;
  double xden;
  const double xinf = 1.79E+308;
  double xm1;
  double xm2;
  double xm4;
  double xnum;
  double y;
  double ysq;

  y = x;

  if ( 0.0 < y && y <= xbig )
  {
    if ( y <= r8_epsilon )
    {
      res = - log ( y );
    }



    else if ( y <= 1.5 )
    {
      if ( y < 0.6796875 )
      {
        corr = -log ( y );
        xm1 = y;
      }
      else
      {
        corr = 0.0;
        xm1 = ( y - 0.5 ) - 0.5;
      }

      if ( y <= 0.5 || 0.6796875 <= y )
      {
        xden = 1.0;
        xnum = 0.0;
        for ( i = 0; i < 8; i++ )
        {
          xnum = xnum * xm1 + p1[i];
          xden = xden * xm1 + q1[i];
        }
        res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) );
      }
      else
      {
        xm2 = ( y - 0.5 ) - 0.5;
        xden = 1.0;
        xnum = 0.0;
        for ( i = 0; i < 8; i++ )
        {
          xnum = xnum * xm2 + p2[i];
          xden = xden * xm2 + q2[i];
        }
        res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) );
      }
    }



    else if ( y <= 4.0 )
    {
      xm2 = y - 2.0;
      xden = 1.0;
      xnum = 0.0;
      for ( i = 0; i < 8; i++ )
      {
        xnum = xnum * xm2 + p2[i];
        xden = xden * xm2 + q2[i];
      }
      res = xm2 * ( d2 + xm2 * ( xnum / xden ) );
    }



    else if ( y <= 12.0 )
    {
      xm4 = y - 4.0;
      xden = -1.0;
      xnum = 0.0;
      for ( i = 0; i < 8; i++ )
      {
        xnum = xnum * xm4 + p4[i];
        xden = xden * xm4 + q4[i];
      }
      res = d4 + xm4 * ( xnum / xden );
    }



    else
    {
      res = 0.0;

      if ( y <= frtbig )
      {
        res = c[6];
        ysq = y * y;
        for ( i = 0; i < 6; i++ )
        {
          res = res / ysq + c[i];
        }
      }
      res = res / y;
      corr = log ( y );
      res = res + sqrtpi - 0.5 * corr;
      res = res + y * ( corr - 1.0 );
    }
  }



  else
  {
    res = xinf;
  }



  return res;
}


double r8_huge ( )




























{
  double value;

  value = 1.79769313486231571E+308;

  return value;
}


double r8_hypot ( double x, double y )

























{
  double a;
  double b;
  double value;

  if ( fabs ( x ) < fabs ( y ) )
  {
    a = fabs ( y );
    b = fabs ( x );
  }
  else
  {
    a = fabs ( x );
    b = fabs ( y );
  }



  if ( a == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = a * sqrt ( 1.0 + ( b / a ) * ( b / a ) );
  }

  return value;
}


bool r8_is_in_01 ( double a )

























{
  bool value;

  value = ( 0.0 <= a && a <= 1.0 );

  return value;
}


bool r8_is_inf ( double r )

























{
  const double r8_huge = 1.79769313486231571E+308;
  bool value;

  if ( r < 0.0 )
  {
    value = ( r < - r8_huge );
  }
  else
  {
    value = ( r8_huge < r );
  }

  return value;
}


bool r8_is_insignificant ( double r, double s )




























{
  double t;
  double tol;
  bool value;

  value = true;

  t = r + s;
  tol = r8_epsilon ( ) * fabs ( r );

  if ( tol < fabs ( r - t ) )
  {
    value = false;
  }
  
  return value;
}


bool r8_is_integer ( double r )

























{
  const int i4_huge = 2147483647;
  bool value;

  if ( ( double ) ( i4_huge ) < r )
  {
    value = false;
  }
  else if ( r < - ( double ) ( i4_huge ) )
  {
    value = false;
  }
  else if ( r == ( double ) ( ( int ) ( r ) ) )
  {
    value = true;
  }
  else
  {
    value = false;
  }
  return value;
}


bool r8_is_nan ( double r )

























{
  bool value;

  value = ( r != r );

  return value;
}


double r8_log_10 ( double x )































{
  double value;

  if ( x == 0.0 )
  {
    value = - r8_big ( );
  }
  else
  {
    value = log10 ( fabs ( x ) );
  }

  return value;
}


double r8_log_2 ( double x )































{
  double value;

  if ( x == 0.0 )
  {
    value = - r8_big ( );
  }
  else
  {
    value = log ( fabs ( x ) ) / log ( 2.0 );
  }

  return value;
}


double r8_log_b ( double x, double b )

































{
  double value;

  if ( b == 0.0 || b == 1.0 || b == -1.0 )
  {
    value = - r8_big ( );
  }
  else if ( fabs ( x ) == 0.0 )
  {
    value = - r8_big ( );
  }
  else
  {
    value = log ( fabs ( x ) ) / log ( fabs ( b ) );
  }

  return value;
}


void r8_mant ( double x, int &s, double &r, int &l )










































{



  if ( x < 0.0 )
  {
    s = -1;
  }
  else
  {
    s = 1;
  }




  if ( x < 0.0 )
  {
    r = -x;
  }
  else
  {
    r = x;
  }

  l = 0;



  if ( x == 0.0 )
  {
    return;
  }

  while ( 2.0 <= r )
  {
    r = r / 2.0;
    l = l + 1;
  }

  while ( r < 1.0 )
  {
    r = r * 2.0;
    l = l - 1;
  }

  return;
}


double r8_max ( double x, double y )





























{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}


double r8_min ( double x, double y )





























{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}


double r8_mod ( double x, double y )













































{
  double value;

  if ( y == 0.0 )
  {
    cerr << "\n";
    cerr << "R8_MOD - Fatal error!\n";
    cerr << "  R8_MOD ( X, Y ) called with Y = " << y << "\n";
    exit ( 1 );
  }

  value = x - ( ( double ) ( ( int ) ( x / y ) ) ) * y;

  if ( x < 0.0 && 0.0 < value )
  {
    value = value - fabs ( y );
  }
  else if ( 0.0 < x && value < 0.0 )
  {
    value = value + fabs ( y );
  }

  return value;
}


double r8_modp ( double x, double y )





















































{
  double value;

  if ( y == 0.0 )
  {
    cerr << "\n";
    cerr << "R8_MODP - Fatal error!\n";
    cerr << "  R8_MODP ( X, Y ) called with Y = " << y << "\n";
    exit ( 1 );
  }

  value = x - ( ( double ) ( ( int ) ( x / y ) ) ) * y;

  if ( value < 0.0 )
  {
    value = value + fabs ( y );
  }

  return value;
}


double r8_mop ( int i )





























{
  double value;

  if ( ( i % 2 ) == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = -1.0;
  }

  return value;
}


int r8_nint ( double x )






































{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( fabs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( fabs ( x ) + 0.5 );
  }

  return value;
}


double r8_normal_01 ( int &seed )
































{
  double r1;
  double r2;
  const double r8_pi = 3.141592653589793;
  double x;

  r1 = r8_uniform_01 ( seed );
  r2 = r8_uniform_01 ( seed );
  x = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * r8_pi * r2 );

  return x;
}


double r8_normal_ab ( double a, double b, int &seed )


































{
  double value;

  value = a + b * r8_normal_01 ( seed );

  return value;
}


double r8_nth_root ( double x, int n )































{
  double e;
  double value;




  if ( x == 0.0 && n == 0 )
  {
    value = 1.0;
    return value;
  }



  if ( x == 0.0 && n < 0 )
  {
    value = NAN;
    return value;
  }



  if ( x < 0.0 && ( n % 2 ) == 0 && 0 < n )
  {
    value = NAN;
    return value;
  }



  if ( n == 0 )
  {
    value = 1.0;
  }



  else if ( n == 1 )
  {
    value = x;
  }



  else if ( n == -1 )
  {
    value = 1.0 / x;
  }
  else
  {
    e = 1.0 / ( double ) ( abs ( n ) );

    if ( 0.0 < x )
    {
      value = pow ( x, e );
    }
    else if ( x == 0.0 )
    {
      value = 0.0;
    }
    else
    {
      value = - pow ( - x, e );
    }

    if ( n < 0 )
    {
      value = 1.0 / value;
    }
  }

  return value;
}


double r8_pi ( )























{
  const double value = 3.141592653589793;

  return value;
}


double r8_pi_sqrt ( )























{
  const double value = 1.7724538509055160273;

  return value;
}


double r8_power ( double r, int p )



























{
  double value;



  if ( p == 0 )
  {
    value = 1.0;
  }




  else if ( r == 0.0 )
  {
    if ( 0 < p )
    {
      value = 0.0;
    }
    else
    {
      value = pow ( r, p );
    }
  }
  else if ( 1 <= p )
  {
    value = pow ( r, p );
  }
  else
  {
    value = pow ( r, p );
  }

  return value;
}


double r8_power_fast ( double r, int p, int &mults )










































{
  int p_mag;
  int p_sign;
  double r2;
  double value;

  mults = 0;



  if ( r == 1.0 )
  {
    value = 1.0;
    return value;
  }

  if ( r == -1.0 )
  {
    if ( ( p % 2 ) == 1 )
    {
      value = -1.0;
    }
    else
    {
      value = 1.0;
    }
    return value;
  }

  if ( r == 0.0 )
  {
    if ( p <= 0 )
    {
      cerr << "\n";
      cerr << "R8_POWER_FAST - Fatal error!\n";
      cerr << "  Base is zero, and exponent is negative.\n";
      exit ( 1 );
    }

    value = 0.0;
    return value;
  }



  if ( p == -1 )
  {
    value = 1.0 / r;
    mults = mults + 1;
    return value;
  }
  else if ( p == 0 )
  {
    value = 1.0;
    return value;
  }
  else if ( p == 1 )
  {
    value = r;
    return value;
  }



  p_mag = abs ( p );
  p_sign = i4_sign ( p );

  value = 1.0;
  r2 = r;

  while ( 0 < p_mag )
  {
    if ( ( p_mag % 2 ) == 1 )
    {
      value = value * r2;
      mults = mults + 1;
    }

    p_mag = p_mag / 2;
    r2 = r2 * r2;
    mults = mults + 1;
  }

  if ( p_sign == -1 )
  {
    value = 1.0 / value;
    mults = mults + 1;
  }

  return value;
}


void r8_print ( double r, string title )

























{
  cout << title << "  "
       << r << "\n";

  return;
}


double r8_radians ( double degrees )

























{
  const double r8_pi = 3.1415926535897932384626434;
  double value;

  value = degrees * r8_pi / 180.0;

  return value;
}


double r8_reverse_bytes ( double x )

























{
  char c;
  union
  {
    double ydouble;
    char ychar[8];
  } y;

  y.ydouble = x;

  c = y.ychar[0];
  y.ychar[0] = y.ychar[7];
  y.ychar[7] = c;

  c = y.ychar[1];
  y.ychar[1] = y.ychar[6];
  y.ychar[6] = c;

  c = y.ychar[2];
  y.ychar[2] = y.ychar[5];
  y.ychar[5] = c;

  c = y.ychar[3];
  y.ychar[3] = y.ychar[4];
  y.ychar[4] = c;

  return ( y.ydouble );
}


double r8_rise ( double x, int n )














































{
  int i;
  double value;

  value = 1.0;

  if ( 0 < n )
  {
    for ( i = 1; i <= n; i++ )
    {
      value = value * x;
      x = x + 1.0;
    }
  }
  else if ( n < 0 )
  {
    for ( i = -1; n <= i; i-- )
    {
      value = value * x;
      x = x - 1.0;
    }
  }

  return value;
}


void r8_rise_values ( int &n_data, double &x, int &n, double &f )























































{
# define N_MAX 15

  static double f_vec[N_MAX] = {
    1680.000000000000,
    1962.597656250000,
    2279.062500000000,
    2631.972656250000,
    3024.000000000000,
    1.000000000000000,
    7.500000000000000,
    63.75000000000000,
    605.6250000000000,
    6359.062500000000,
    73129.21875000000,
    914115.2343750000,
    1.234055566406250E+07,
    1.789380571289063E+08,
    2.773539885498047E+09 };

  static int n_vec[N_MAX] = {
    4,
    4,
    4,
    4,
    4,
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9 };

  static double x_vec[N_MAX] = {
    5.00,
    5.25,
    5.50,
    5.75,
    6.00,
    7.50,
    7.50,
    7.50,
    7.50,
    7.50,
    7.50,
    7.50,
    7.50,
    7.50,
    7.50 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    n = 0;
    f = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    n = n_vec[n_data-1];
    f = f_vec[n_data-1];
  }

  return;
# undef N_MAX
}


double r8_round ( double x )






































{
  double value;

  if ( x < 0.0 )
  {
    value = - ( double ) floor ( - x + 0.5 );
  }
  else
  {
    value =   ( double ) floor (   x + 0.5 );
  }

  return value;
}


int r8_round_i4 ( double x )






































{
  int value;

  if ( x < 0.0 )
  {
    value = - floor ( - x + 0.5 );
  }
  else
  {
    value =   floor (   x + 0.5 );
  }

  return value;
}


double r8_round2 ( int nplace, double x )





















































{
  int iplace;
  int l;
  int s;
  double xmant;
  double xtemp;
  double value;

  value = 0.0;



  if ( x == 0.0 )
  {
    return value;
  }

  if ( nplace <= 0 )
  {
    return value;
  }



  if ( 0.0 < x )
  {
    s = 1;
    xtemp = x;
  }
  else
  {
    s = -1;
    xtemp = -x;
  }




  l = 0;

  while ( 2.0 <= xtemp )
  {
    xtemp = xtemp / 2.0;
    l = l + 1;
  }

  while ( xtemp < 1.0 )
  {
    xtemp = xtemp * 2.0;
    l = l - 1;
  }



  xmant = 0.0;
  iplace = 0;

  for ( ; ; )
  {
    xmant = 2.0 * xmant;

    if ( 1.0 <= xtemp )
    {
      xmant = xmant + 1.0;
      xtemp = xtemp - 1.0;
    }

    iplace = iplace + 1;

    if ( xtemp == 0.0 || nplace <= iplace )
    {
      value = s * xmant * pow ( 2.0, l );
      break;
    }

    l = l - 1;
    xtemp = xtemp * 2.0;
  }

  return value;
}


double r8_roundb ( int base, int nplace, double x )

































































{
  int iplace;
  int is;
  int js;
  int l;
  double r8_base;
  double value;
  double xmant;
  double xtemp;

  value = 0.0;
  r8_base = ( double ) base;



  if ( base == 0 )
  {
    cerr << "\n";
    cerr << "R8_ROUNDB - Fatal error!\n";
    cerr << "  The base BASE cannot be zero.\n";
    exit ( 1 );
  }



  if ( x == 0.0 )
  {
    return value;
  }

  if ( nplace <= 0 )
  {
    return value;
  }



  if ( 0.0 < x )
  {
    is = 1;
    xtemp = x;
  }
  else
  {
    is = -1;
    xtemp = -x;
  }




  l = 0;

  while ( fabs ( r8_base ) <= fabs ( xtemp ) )
  {
    xtemp = xtemp / r8_base;

    if ( xtemp < 0.0 )
    {
      is = -is;
      xtemp = -xtemp;
    }
    l = l + 1;
  }

  while ( fabs ( xtemp ) < 1.0 )
  {
    xtemp = xtemp * r8_base;

    if ( xtemp < 0.0 )
    {
      is = -is;
      xtemp = -xtemp;
    }

    l = l - 1;
  }




  xmant = 0.0;
  iplace = 0;
  js = is;

  for ( ; ; )
  {
    xmant = r8_base * xmant;

    if ( xmant < 0.0 )
    {
      js = -js;
      xmant = -xmant;
    }

    if ( 1.0 <= xtemp )
    {
      xmant = xmant + ( int ) ( xtemp );
      xtemp = xtemp - ( int ) ( xtemp );
    }

    iplace = iplace + 1;

    if ( xtemp == 0.0 || nplace <= iplace )
    {
      value = ( double ) js * xmant * pow ( r8_base, l );
      break;
    }

    l = l - 1;
    xtemp = xtemp * r8_base;

    if ( xtemp < 0.0 )
    {
      is = -is;
      xtemp = -xtemp;
    }
  }

  return value;
}


double r8_roundx ( int nplace, double x )


























































{
  int iplace;
  int is;
  int l;
  double xmant;
  double xround;
  double xtemp;

  xround = 0.0;



  if ( x == 0.0 )
  {
    return xround;
  }

  if ( nplace <= 0 )
  {
    return xround;
  }



  if ( 0.0 < x )
  {
    is = 1;
    xtemp = x;
  }
  else
  {
    is = -1;
    xtemp = -x;
  }




  l = 0;

  while ( 10.0 <= x )
  {
    xtemp = xtemp / 10.0;
    l = l + 1;
  }

  while ( xtemp < 1.0 )
  {
    xtemp = xtemp * 10.0;
    l = l - 1;
  }




  xmant = 0.0;
  iplace = 0;

  for ( ; ; )
  {
    xmant = 10.0 * xmant;

    if ( 1.0 <= xtemp )
    {
      xmant = xmant + ( int ) xtemp;
      xtemp = xtemp - ( int ) xtemp;
    }

    iplace = iplace + 1;

    if ( xtemp == 0.0 || nplace <= iplace )
    {
      xround = is * xmant * pow ( 10.0, l );
      break;
    }

    l = l - 1;
    xtemp = xtemp * 10.0;
  }

  return xround;
}


double r8_secd ( double degrees )

























{
  const double r8_pi = 3.141592653589793;
  double radians;
  double value;

  radians = r8_pi * ( degrees / 180.0 );

  value = 1.0 / cos ( radians );

  return value;
}


double r8_sech ( double x )

























{
  const double log_huge = 80.0;
  double value;

  if ( log_huge < fabs ( x ) )
  {
    value = 0.0;
  }
  else
  {
    value = 1.0 / cosh ( x );
  }
  return value;
}


double r8_sign ( double x )

























{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}


double r8_sign3 ( double x )

























{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else if ( x == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}


char r8_sign_char ( double x )

























{
  char value;

  if ( x < 0.0 )
  {
    value = '-';
  }
  else if ( x == 0.0 )
  {
    value = '0';
  }
  else
  {
    value = '+';
  }
  return value;
}


bool r8_sign_match ( bool r1, bool r2 )
































{
  bool value;

  value = ( r1 <= 0.0 && r2 <= 0.0 ) || ( 0.0 <= r1 && 0.0 <= r2 );

  return value;
}


bool r8_sign_match_strict ( bool r1, bool r2 )

























{
  bool value;

  value = ( r1 < 0.0 && r2 < 0.0 ) || 
          ( r1 == 0.0 && r2 == 0.0 ) || 
          ( 0.0 < r1 && 0.0 < r2 );

  return value;
}


bool r8_sign_opposite ( double r1, double r2 )
































{
  bool value;

  value = ( r1 <= 0.0 && 0.0 <= r2 ) || ( r2 <= 0.0 && 0.0 <= r1 );

  return value;
}


bool r8_sign_opposite_strict ( double r1, double r2 )
































{
  bool value;

  value = ( r1 < 0.0 && 0.0 < r2 ) || ( r2 < 0.0 && 0.0 < r1 );

  return value;
}


double r8_sign2 ( double x, double y )


























{
  double value;

  if ( 0.0 <= y )
  {
    value = fabs ( x );
  } 
  else
  {
    value = - fabs ( x );
  }
  return value;
}


void r8_sincos_sum ( double a, double b, double &d, double &e, double &f )



































{
  const double r8_pi = 3.141592653589793E+00;

  d = sqrt ( a * a + b * b );
  e = atan2 ( b, a );
  f = atan2 ( b, a ) - r8_pi / 2.0E+00;
  if ( f < - r8_pi )
  {
    f = f + 2.0E+00 * r8_pi;
  }

  return;
}


double r8_sind ( double degrees )

























{
  const double r8_pi = 3.141592653589793;
  double radians;
  double value;

  radians = r8_pi * ( degrees / 180.0 );

  value = sin ( radians );

  return value;
}


double r8_sqrt_i4 ( int i )

























{
  double value;

  value = sqrt ( ( double ) ( i ) );

  return value;
}


double r8_sum ( double x, double y )

























{
  double value;

  value = x + y;

  return value;
}


void r8_swap ( double &x, double &y )
























{
  double z;

  z = x;
  x = y;
  y = z;

  return;
}


void r8_swap3 ( double &x, double &y, double &z )

































{
  double w;

  w = x;
  x = y;
  y = z;
  z =  w;

  return;
}


double r8_tand ( double degrees )

























{
  const double r8_pi = 3.141592653589793;
  double radians;
  double value;

  radians = r8_pi * ( degrees / 180.0 );

  value = sin ( radians ) / cos ( radians );

  return value;
}


double r8_tiny ( )























{
  const double value = 0.4450147717014E-307;

  return value;
}


void r8_to_dhms ( double r, int &d, int &h, int &m, int &s )


























{
  int sign;

  if ( 0.0 <= r )
  {
    sign = 1;
  }
  else
  {
    sign = - 1;
    r = - r;
  }

  d = ( int ) r;

  r = r - ( double ) d;
  r = 24.0 * r;
  h = ( int ) r;

  r = r - ( double ) h;
  r = 60.0 * r;
  m = ( int ) r;

  r = r - ( double ) m;
  r = 60.0 * r;
  s = ( int ) r;

  if ( sign == -1 )
  {
    d = -d;
    h = -h;
    m = -m;
    s = -s;
  }

  return;
}


int r8_to_i4 ( double xmin, double xmax, double x, int ixmin, int ixmax )







































{
  int ix;
  double temp;

  if ( xmax == xmin )
  {
    cerr << "\n";
    cerr << "R8_TO_I4 - Fatal error!\n";
    cerr << "  XMAX = XMIN, making a zero divisor.\n";
    cerr << "  XMAX = " << xmax << "\n";
    cerr << "  XMIN = " << xmin << "\n";
    exit ( 1 );
  }

  temp =
      ( ( xmax - x        ) * ( double ) ixmin
      + (        x - xmin ) * ( double ) ixmax )
      / ( xmax     - xmin );

  if ( 0.0 <= temp )
  {
    temp = temp + 0.5;
  }
  else
  {
    temp = temp - 0.5;
  }

  ix = ( int ) temp;

  return ix;
}


double r8_to_r8_discrete ( double r, double rmin, double rmax, int nr )













































{
  int f;
  double rd;



  if ( nr < 1 )
  {
    cerr << "\n";
    cerr << "R8_TO_R8_DISCRETE - Fatal error!\n";
    cerr << "  NR = " << nr << "\n";
    cerr << "  but NR must be at least 1.\n";
    exit ( 1 );
  }

  if ( nr == 1 )
  {
    rd = 0.5 * ( rmin + rmax );
    return rd;
  }

  if ( rmax == rmin )
  {
    rd = rmax;
    return rd;
  }

  f = r8_nint ( ( double ) ( nr ) * ( rmax - r ) / ( rmax - rmin ) );
  f = i4_max ( f, 0 );
  f = i4_min ( f, nr );

  rd = ( ( double ) (      f ) * rmin
       + ( double ) ( nr - f ) * rmax )
       / ( double ) ( nr     );

  return rd;
}


double r8_uniform_01 ( int &seed )











































































{
  const int i4_huge = 2147483647;
  int k;
  double r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }
  r = ( double ) ( seed ) * 4.656612875E-10;

  return r;
}


double r8_uniform_ab ( double a, double b, int &seed )

































{
  const int i4_huge = 2147483647;
  int k;
  double value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  value = ( double ) ( seed ) * 4.656612875E-10;

  value = a + ( b - a ) * value;

  return value;
}


void r8_unswap3 ( double &x, double &y, double &z )

































{
  double w;

  w = z;
  z = y;
  y = x;
  x = w;

  return;
}


double r8_walsh_1d ( double x, int digit )



































{
  int n;
  double value;



  x = fabs ( x );




  x = x / pow ( 2.0, digit );





  n = ( int ) x;



  if ( ( n % 2 ) == 0 )
  {
    value = 0.0;
  }
  else
  {
    value = 1.0;
  }

  return value;
}


double r8_wrap ( double r, double rlo, double rhi )























































{
  int n;
  double rhi2;
  double rlo2;
  double rwide;
  double value;



  if ( rlo <= rhi )
  {
    rlo2 = rlo;
    rhi2 = rhi;
  }
  else
  {
    rlo2 = rhi;
    rhi2 = rlo;
  }



  rwide = rhi2 - rlo2;




  if ( rwide == 0.0 )
  {
    value = rlo;
  }
  else if ( r < rlo2 )
  {
    n = ( int ) ( ( rlo2 - r ) / rwide ) + 1;
    value = r + n * rwide;
    if ( value == rhi )
    {
      value = rlo;
    }
  }
  else
  {
    n = ( int ) ( ( r - rlo2 ) / rwide );
    value = r - n * rwide;
    if ( value == rlo )
    {
      value = rhi;
    }
  }
  return value;
}


double r82_dist_l2 ( double a1[2], double a2[2] )

































{
  double value;

  value = sqrt ( pow ( a1[0] - a2[0], 2 )
               + pow ( a1[1] - a2[1], 2 ) );

  return value;
}


void r82_print ( double a[2], string title )



































{
  cout << "  " << title << " : ";
  cout << ": ( " << setw(12) << a[0]
       << ", "   << setw(12) << a[1] << " )\n";

  return;
}


void r82_uniform_ab ( double b, double c, int &seed, double r[] )































{
  int i;

  for ( i = 0; i < 2; i++ )
  {
    r[i] = r8_uniform_ab ( b, c, seed );
  }

  return;
}


void r82col_print_part ( int n, double a[], int max_print, string title )











































{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(14) << a[i+0*n]
           << "  " << setw(14) << a[i+1*n] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[i+0*n]
           << "  " << setw(14) << a[i+1*n]  << "\n";
    }
    cout << "  ........  ..............  ..............\n";
    i = n - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[i+0*n]
         << "  " << setw(14) << a[i+1*n]  << "\n";
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[i+0*n]
           << "  " << setw(14) << a[i+1*n]  << "\n";
    }
    i = max_print - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[i+0*n]
         << "  " << setw(14) << a[i+1*n] 
         << "  " << "...more entries...\n";
  }

  return;
}


void r82poly2_print ( double a, double b, double c, double d, double e,
  double f )























{
  cout << "  " << setw(8) << a
       << " * x^2 + " << setw(8) << b
       << " * y^2 + " << setw(8) << c
       << " * xy  + " << "\n";
  cout << "  " << setw(8) << d
       << " * x   + " << setw(8) << e
       << " * y   + " << setw(8) << f << "\n";

  return;
}


int r82poly2_type ( double a, double b, double c, double d, double e, double f )
































































{
  double delta;
  double j;
  double k;
  int type = 0;



  if ( a == 0.0 && b == 0.0 && c == 0.0 )
  {
    if ( d == 0.0 && e == 0.0 )
    {
      if ( f == 0.0 )
      {
        type = 11;
      }
      else
      {
        type = 12;
      }
    }
    else
    {
      type = 10;
    }
    return type;
  }

  delta =
      8.0 * a * b * f
    + 2.0 * c * e * d
    - 2.0 * a * e * e
    - 2.0 * b * d * d
    - 2.0 * f * c * c;

  j = 4.0 * a * b - c * c;

  if ( delta != 0.0 )
  {
    if ( j < 0.0 )
    {
      type = 1;
    }
    else if ( j == 0.0 )
    {
      type = 2;
    }
    else if ( 0.0 < j )
    {
      if ( r8_sign ( delta ) != r8_sign ( a + b ) )
      {
        type = 3;
      }
      else if ( r8_sign ( delta ) == r8_sign ( a + b ) )
      {
        type = 4;
      }
    }
  }
  else if ( delta == 0.0 )
  {
    if ( j < 0.0 )
    {
      type = 5;
    }
    else if ( 0.0 < j )
    {
      type = 6;
    }
    else if ( j == 0.0 )
    {
      k = 4.0 * ( a + b ) * f - d * d - e * e;

      if ( k < 0.0 )
      {
        type = 7;
      }
      else if ( 0.0 < k )
      {
        type = 8;
      }
      else if ( k == 0.0 )
      {
        type = 9;
      }
    }
  }

  return type;
}


void r82poly2_type_print ( int type )























{
  if ( type == 1 )
  {
    cout << "  The set of solutions forms a hyperbola.\n";
  }
  else if ( type == 2 )
  {
    cout << "  The set of solutions forms a parabola.\n";
  }
  else if ( type == 3 )
  {
    cout << "  The set of solutions forms an ellipse.\n";
  }
  else if ( type == 4 )
  {
    cout << "  The set of solutions forms an imaginary ellipse.\n";
    cout << "  (There are no real solutions).\n";
  }
  else if ( type == 5 )
  {
    cout << "  The set of solutions forms a pair of intersecting lines.\n";
  }
  else if ( type == 6 )
  {
    cout << "  The set of solutions is a single point.\n";
  }
  else if ( type == 7 )
  {
    cout << "  The set of solutions form a pair of distinct parallel lines.\n";
  }
  else if ( type == 8 )
  {
    cout << "  The set of solutions forms a pair of imaginary parallel lines.\n";
    cout << "  (There are no real solutions).\n";
  }
  else if ( type == 9 )
  {
    cout << "  The set of solutions forms a pair of coincident lines.\n";
  }
  else if ( type == 10 )
  {
    cout << "  The set of solutions forms a single line.\n";
  }
  else if ( type == 11 )
  {
    cout << "  The set of solutions is all space.\n";
  }
  else if ( type == 12 )
  {
    cout << "  The set of solutions is empty.\n";
  }
  else
  {
    cout << "  This type index is unknown.\n";
  }
  return;
}


double *r82row_max ( int n, double a[] )































{
# define DIM_NUM 2

  double *amax = NULL;
  int i;
  int j;

  if ( n <= 0 )
  {
    return NULL;
  }

  amax = new double[DIM_NUM];

  for ( i = 0; i < DIM_NUM; i++ )
  {
    amax[i] = a[i+0*DIM_NUM];
    for ( j = 1; j < n; j++ )
    {
      if ( amax[i] < a[0+j*DIM_NUM] )
      {
        amax[i] = a[0+j*DIM_NUM];
      }
    }
  }
  return amax;
# undef DIM_NUM
}


double *r82row_min ( int n, double a[] )































{
# define DIM_NUM 2

  double *amin = NULL;
  int i;
  int j;

  if ( n <= 0 )
  {
    return NULL;
  }

  amin = new double[DIM_NUM];

  for ( i = 0; i < DIM_NUM; i++ )
  {
    amin[i] = a[i+0*DIM_NUM];
    for ( j = 1; j < n; j++ )
    {
      if ( a[0+j*DIM_NUM] < amin[i] )
      {
        amin[i] = a[0+j*DIM_NUM];
      }
    }
  }
  return amin;
# undef DIM_NUM
}


int r82row_order_type ( int n, double a[] )









































{
  int i;
  int order;



  i = 0;

  for ( ; ; )
  {
    i = i + 1;

    if ( n <= i )
    {
      order = 0;
      return order;
    }

    if ( a[0+0*2] < a[0+i*2] || ( a[0+0*2] == a[0+i*2] && a[1+0*2] < a[1+i*2] ) )
    {
      if ( i == 2 )
      {
        order = 2;
      }
      else
      {
        order = 1;
      }
      break;
    }
    else if ( a[0+i*2] < a[0+0*2] || 
      ( a[0+i*2] == a[0+0*2] && a[1+i*2] < a[1+0*2] ) )
    {
      if ( i == 2 )
      {
        order = 4;
      }
      else
      {
        order = 3;
      }
      break;
    }
  }



  for ( ; ; )
  {
    i = i + 1;
    if ( n <= i )
    {
      break;
    }

    if ( order == 1 )
    {
      if ( a[0+i*2] < a[0+(i-1)*2] ||
        ( a[0+i*2] == a[0+(i-1)*2] && a[1+i*2] < a[1+(i-1)*2] ) )
      {
        order = -1;
        break;
      }
    }
    else if ( order == 2 )
    {
      if ( a[0+i*2] < a[0+(i-1)*2] ||
        ( a[0+i*2] == a[0+(i-1)*2] && a[1+i*2] < a[1+(i-1)*2] ) )
      {
        order = -1;
        break;
      }
      else if ( a[0+i*2] == a[0+(i-1)*2] && a[1+i*2] == a[1+(i-1)*2] )
      {
        order = 1;
      }
    }
    else if ( order == 3 )
    {
      if ( a[0+(i-1)*2] < a[0+i*2] ||
        ( a[0+(i-1)*2] == a[0+i*2] && a[1+(i-1)*2] < a[1+i*2] ) )
      {
        order = -1;
        break;
      }
    }
    else if ( order == 4 )
    {
      if ( a[0+(i-1)*2] < a[0+i*2] ||
        ( a[0+(i-1)*2] == a[0+i*2] && a[1+(i-1)*2] < a[1+i*2] ) )
      {
        order = -1;
        break;
      }
      else if ( a[0+i*2] == a[0+(i-1)*2] && a[1+i*2] == a[1+(i-1)*2] )
      {
        order = 3;
      }
    }
  }
  return order;
}


void r82row_part_quick_a ( int n, double a[], int &l, int &r )
























































{
  int i;
  int j;
  double key[2];
  int ll;
  int m;
  int rr;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R82ROW_PART_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    l = 0;
    r = 2;
    return;
  }

  key[0] = a[2*0+0];
  key[1] = a[2*0+1];
  m = 1;



  ll = 1;
  rr = n + 1;

  for ( i = 2; i <= n; i++ )
  {
    if ( r8vec_gt ( 2, a+2*ll, key ) )
    {
      rr = rr - 1;
      r8vec_swap ( 2, a+2*(rr-1), a+2*ll );
    }
    else if ( r8vec_eq ( 2, a+2*ll, key ) )
    {
      m = m + 1;
      r8vec_swap ( 2, a+2*(m-1), a+2*ll );
      ll = ll + 1;
    }
    else if ( r8vec_lt ( 2, a+2*ll, key ) )
    {
      ll = ll + 1;
    }

  }



  for ( i = 0; i < ll - m; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      a[2*i+j] = a[2*(i+m)+j];
    }
  }

  ll = ll - m;

  for ( i = ll; i < ll+m; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      a[2*i+j] = key[j];
    }
  }

  l = ll;
  r = rr;

  return;
}


void r82row_permute ( int n, int p[], double a[] )






















































{
  double a_temp[2];
  int i;
  int iget;
  int iput;
  int istart;

  if ( !perm0_check ( n, p ) )
  {
    cerr << "\n";
    cerr << "R82ROW_PERMUTE - Fatal error!\n";
    cerr << "  PERM0_CHECK rejects permutation.\n";
    exit ( 1 );
  }





  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1;
  }



  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      a_temp[0] = a[0+(istart-1)*2];
      a_temp[1] = a[1+(istart-1)*2];
      iget = istart;



      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cerr << "\n";
          cerr << "R82ROW_PERMUTE - Fatal error!\n";
          cerr << "  Entry IPUT = " << iput << " of the permutation has\n";
          cerr << "  an illegal value IGET = " << iget << ".\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[0+(iput-1)*2] = a_temp[0];
          a[1+(iput-1)*2] = a_temp[1];
          break;
        }
        a[0+(iput-1)*2] = a[0+(iget-1)*2];
        a[1+(iput-1)*2] = a[1+(iget-1)*2];
      }
    }
  }



  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }



  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1;
  }
  return;
}


void r82row_print ( int n, double a[], string title )































{
  int j;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    cout << "  " << setw(8)  << j
         << ": " << setw(14) << a[0+j*2]
         << "  " << setw(14) << a[1+j*2] << "\n";
  }

  return;
}


void r82row_print_part ( int n, double a[], int max_print, string title )











































{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(14) << a[0+i*2]
           << "  " << setw(14) << a[1+i*2] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[0+i*2]
           << "  " << setw(14) << a[1+i*2]  << "\n";
    }
    cout << "  ........  ..............  ..............\n";
    i = n - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[0+i*2]
         << "  " << setw(14) << a[1+i*2]  << "\n";
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[0+i*2]
           << "  " << setw(14) << a[1+i*2]  << "\n";
    }
    i = max_print - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[0+i*2]
         << "  " << setw(14) << a[1+i*2] 
         << "  " << "...more entries...\n";
  }

  return;
}


int *r82row_sort_heap_index_a ( int n, double a[] )
















































{
  double aval[2];
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = new int[n];

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    indx[0] = indx[0];
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval[0] = a[0+indxt*2];
      aval[1] = a[1+indxt*2];
    }
    else
    {
      indxt = indx[ir-1];
      aval[0] = a[0+indxt*2];
      aval[1] = a[1+indxt*2];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }
    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if (   a[0+indx[j-1]*2] <  a[0+indx[j]*2] ||
             ( a[0+indx[j-1]*2] == a[0+indx[j]*2] &&
               a[1+indx[j-1]*2] <  a[1+indx[j]*2] ) )
        {
          j = j + 1;
        }
      }

      if (   aval[0] <  a[0+indx[j-1]*2] ||
           ( aval[0] == a[0+indx[j-1]*2] &&
             aval[1] <  a[1+indx[j-1]*2] ) )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return indx;
}


void r82row_sort_quick_a ( int n, double a[] )































{
# define LEVEL_MAX 30

  int base;
  int l_segment;
  int level;
  int n_segment;
  int rsave[LEVEL_MAX];
  int r_segment;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R82ROW_SORT_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    return;
  }

  level = 1;
  rsave[level-1] = n + 1;
  base = 1;
  n_segment = n;

  while ( 0 < n_segment )
  {



    r82row_part_quick_a ( n_segment, a+2*(base-1)+0, l_segment, r_segment );



    if ( 1 < l_segment )
    {
      if ( LEVEL_MAX < level )
      {
        cerr << "\n";
        cerr<< "R82ROW_SORT_QUICK_A - Fatal error!\n";
        cerr << "  Exceeding recursion maximum of " << LEVEL_MAX << "\n";
        exit ( 1 );
      }

      level = level + 1;
      n_segment = l_segment;
      rsave[level-1] = r_segment + base - 1;
    }




    else if ( r_segment < n_segment )
    {
      n_segment = n_segment + 1 - r_segment;
      base = base + r_segment - 1;
    }



    else
    {
      for ( ; ; )
      {
        if ( level <= 1 )
        {
          n_segment = 0;
          break;
        }

        base = rsave[level-1];
        n_segment = rsave[level-2] - rsave[level-1];
        level = level - 1;

        if ( 0 < n_segment )
        {
          break;
        }
      }
    }
  }
  return;
# undef LEVEL_MAX
}


double r83_norm ( double x, double y, double z )





























{
  double value;

  value = sqrt ( x * x + y * y + z * z );

  return value;
}


void r83col_print_part ( int n, double a[], int max_print, string title )











































{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(14) << a[i+0*n]
           << "  " << setw(14) << a[i+1*n] 
           << "  " << setw(14) << a[i+2*n] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[i+0*n]
           << "  " << setw(14) << a[i+1*n] 
           << "  " << setw(14) << a[i+2*n]  << "\n";
    }
    cout << "  ........  ..............  ..............  ..............\n";
    i = n - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[i+0*n]
         << "  " << setw(14) << a[i+1*n] 
         << "  " << setw(14) << a[i+2*n]  << "\n";
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[i+0*n]
           << "  " << setw(14) << a[i+1*n] 
           << "  " << setw(14) << a[i+2*n]  << "\n";
    }
    i = max_print - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[i+0*n]
         << "  " << setw(14) << a[i+1*n] 
         << "  " << setw(14) << a[i+2*n] 
         << "  " << "...more entries...\n";
  }

  return;
}


double *r83row_max ( int n, double a[] )































{
# define DIM_NUM 3

  double *amax = NULL;
  int i;
  int j;

  if ( n <= 0 )
  {
    return NULL;
  }

  amax = new double[DIM_NUM];

  for ( i = 0; i < DIM_NUM; i++ )
  {
    amax[i] = a[i+0*DIM_NUM];
    for ( j = 1; j < n; j++ )
    {
      if ( amax[i] < a[i+j*DIM_NUM] )
      {
        amax[i] = a[i+j*DIM_NUM];
      }
    }
  }
  return amax;
# undef DIM_NUM
}


double *r83row_min ( int n, double a[] )































{
# define DIM_NUM 3

  double *amin = NULL;
  int i;
  int j;

  if ( n <= 0 )
  {
    return NULL;
  }

  amin = new double[DIM_NUM];

  for ( i = 0; i < DIM_NUM; i++ )
  {
    amin[i] = a[i+0*DIM_NUM];
    for ( j = 1; j < n; j++ )
    {
      if ( a[i+j*DIM_NUM] < amin[i] )
      {
        amin[i] = a[i+j*DIM_NUM];
      }
    }
  }
  return amin;
# undef DIM_NUM
}


void r83row_part_quick_a ( int n, double a[], int &l, int &r )








































{
  int i;
  int j;
  double key[3];
  int ll;
  int m;
  int rr;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R83ROW_PART_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    l = 0;
    r = 2;
    return;
  }

  key[0] = a[3*0+0];
  key[1] = a[3*0+1];
  key[2] = a[3*0+2];
  m = 1;



  ll = 1;
  rr = n + 1;

  for ( i = 2; i <= n; i++ )
  {
    if ( r8vec_gt ( 3, a+3*ll, key ) )
    {
      rr = rr - 1;
      r8vec_swap ( 3, a+3*(rr-1), a+3*ll );
    }
    else if ( r8vec_eq ( 3, a+3*ll, key ) )
    {
      m = m + 1;
      r8vec_swap ( 3, a+3*(m-1), a+3*ll );
      ll = ll + 1;
    }
    else if ( r8vec_lt ( 3, a+3*ll, key ) )
    {
      ll = ll + 1;
    }
  }



  for ( i = 0; i < ll - m; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      a[3*i+j] = a[3*(i+m)+j];
    }
  }

  ll = ll - m;

  for ( i = ll; i < ll+m; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      a[3*i+j] = key[j];
    }
  }

  l = ll;
  r = rr;

  return;
}


void r83row_print_part ( int n, double a[], int max_print, string title )











































{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(14) << a[0+i*3]
           << "  " << setw(14) << a[1+i*3] 
           << "  " << setw(14) << a[2+i*3] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[0+i*3]
           << "  " << setw(14) << a[1+i*3] 
           << "  " << setw(14) << a[2+i*3]  << "\n";
    }
    cout << "  ........  ..............  ..............  ..............\n";
    i = n - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[0+i*3]
         << "  " << setw(14) << a[1+i*3] 
         << "  " << setw(14) << a[2+i*3]  << "\n";
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[0+i*3]
           << "  " << setw(14) << a[1+i*3] 
           << "  " << setw(14) << a[2+i*3]  << "\n";
    }
    i = max_print - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[0+i*3]
         << "  " << setw(14) << a[1+i*3] 
         << "  " << setw(14) << a[2+i*3] 
         << "  " << "...more entries...\n";
  }

  return;
}


void r83row_sort_quick_a ( int n, double a[] )































{
# define LEVEL_MAX 30

  int base;
  int l_segment;
  int level;
  int n_segment;
  int rsave[LEVEL_MAX];
  int r_segment;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R83ROW_SORT_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    return;
  }

  level = 1;
  rsave[level-1] = n + 1;
  base = 1;
  n_segment = n;

  while ( 0 < n_segment )
  {



    r83row_part_quick_a ( n_segment, a+3*(base-1)+0, l_segment, r_segment );



    if ( 1 < l_segment )
    {
      if ( LEVEL_MAX < level )
      {
        cerr << "\n";
        cerr << "R83ROW_SORT_QUICK_A - Fatal error!\n";
        cerr << "  Exceeding recursion maximum of " << LEVEL_MAX << "\n";
        exit ( 1 );
      }
      level = level + 1;
      n_segment = l_segment;
      rsave[level-1] = r_segment + base - 1;
    }




    else if ( r_segment < n_segment )
    {
      n_segment = n_segment + 1 - r_segment;
      base = base + r_segment - 1;
    }



    else
    {
      for ( ; ; )
      {
        if ( level <= 1 )
        {
          n_segment = 0;
          break;
        }

        base = rsave[level-1];
        n_segment = rsave[level-2] - rsave[level-1];
        level = level - 1;

        if ( 0 < n_segment )
        {
          break;
        }
      }
    }
  }
  return;
# undef LEVEL_MAX
}


void r8block_delete ( int l, int m, int n, double ***a )
































{
  int i;
  int j;

  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      delete [] a[i][j];
    }
  }

  for ( i = 0; i < l; i++ )
  {
    delete [] a[i];
  }

  delete [] a;

  return;
}


double *r8block_expand_linear ( int l, int m, int n, double x[], int lfat,
  int mfat, int nfat )











































{
  int i;
  int ihi;
  int ii;
  int iii;
  int ip1;
  int j;
  int jhi;
  int jj;
  int jjj;
  int jp1;
  int k;
  int khi;
  int kk;
  int kkk;
  int kp1;
  int l2;
  int m2;
  int n2;
  double r;
  double s;
  double t;
  double x000;
  double x001;
  double x010;
  double x011;
  double x100;
  double x101;
  double x110;
  double x111;
  double *xfat;

  l2 = ( l - 1 ) * ( lfat + 1 ) + 1;
  m2 = ( m - 1 ) * ( mfat + 1 ) + 1;
  n2 = ( n - 1 ) * ( nfat + 1 ) + 1;

  xfat = new double[l2*m2*n2];

  for ( i = 1; i <= l; i++ )
  {
    if ( i < l )
    {
      ihi = lfat;
    }
    else
    {
      ihi = 0;
    }

    for ( j = 1; j <= m; j++ )
    {
      if ( j < m )
      {
        jhi = mfat;
      }
      else
      {
        jhi = 0;
      }

      for ( k = 1; k <= n; k++ )
      {
        if ( k < n )
        {
          khi = nfat;
        }
        else
        {
          khi = 0;
        }

        if ( i < l )
        {
          ip1 = i + 1;
        }
        else
        {
          ip1 = i;
        }

        if ( j < m )
        {
          jp1 = j + 1;
        }
        else
        {
          jp1 = j;
        }

        if ( k < n )
        {
          kp1 = k + 1;
        }
        else
        {
          kp1 = k;
        }

        x000 = x[i-1+(j-1)*l+(k-1)*l*m];
        x001 = x[i-1+(j-1)*l+(kp1-1)*l*m];
        x100 = x[ip1-1+(j-1)*l+(k-1)*l*m];
        x101 = x[ip1-1+(j-1)*l+(kp1-1)*l*m];
        x010 = x[i-1+(jp1-1)*l+(k-1)*l*m];
        x011 = x[i-1+(jp1-1)*l+(kp1-1)*l*m];
        x110 = x[ip1-1+(jp1-1)*l+(k-1)*l*m];
        x111 = x[ip1-1+(jp1-1)*l+(kp1-1)*l*m];

        for ( ii = 0; ii <= ihi; ii++ )
        {
          r = ( double ) ( ii ) / ( double ) ( ihi + 1 );

          for ( jj = 0; jj <= jhi; jj++ )
          {
            s = ( double ) ( jj ) / ( double ) ( jhi + 1 );

            for ( kk = 0; kk <= khi; kk++ )
            {
              t = ( double ) ( kk ) / ( double ) ( khi + 1 );

              iii = 1 + ( i - 1 ) * ( lfat + 1 ) + ii;
              jjj = 1 + ( j - 1 ) * ( mfat + 1 ) + jj;
              kkk = 1 + ( k - 1 ) * ( nfat + 1 ) + kk;

              xfat[iii-1+(jjj-1)*l2+(kkk-1)*l2*m2] =
                  x000 * ( 1.0 - r ) * ( 1.0 - s ) * ( 1.0 - t )
                + x001 * ( 1.0 - r ) * ( 1.0 - s ) * (       t )
                + x010 * ( 1.0 - r ) * (       s ) * ( 1.0 - t )
                + x011 * ( 1.0 - r ) * (       s ) * (       t )
                + x100 * (       r ) * ( 1.0 - s ) * ( 1.0 - t )
                + x101 * (       r ) * ( 1.0 - s ) * (       t )
                + x110 * (       r ) * (       s ) * ( 1.0 - t )
                + x111 * (       r ) * (       s ) * (       t );
            }
          }
        }
      }
    }
  }

  return xfat;
}


double ***r8block_new ( int l, int m, int n )





































{
  double ***a;
  int i;
  int j;

  a = new double **[l];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "R8BLOCK_NEW - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( i = 0; i < l; i++ )
  {
    a[i] = new double *[m];
    if ( a[i] == NULL )
    {
      cerr << "\n";
      cerr << "R8BLOCK_NEW - Fatal error!\n";
      cerr << "  Unable to allocate column pointer array.\n";
      exit ( 1 );
    }
  }

  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      a[i][j] = new double[n];
      if ( a[i][j] == NULL )
      {
        cerr << "\n";
        cerr << "R8BLOCK_NEW - Fatal error!\n";
        cerr << "  Unable to allocate layer array.\n";
        exit ( 1 );
      }
    }
  }
  return a;
}


void r8block_print ( int l, int m, int n, double a[], string title )



























{
  int i;
  int j;
  int jhi;
  int jlo;
  int k;

  cout << "\n";
  cout << title << "\n";

  for ( k = 1; k <= n; k++ )
  {
    cout << "\n";
    cout << "  K = " << k << "\n";
    cout << "\n";
    for ( jlo = 1; jlo <= m; jlo = jlo + 5 )
    {
      jhi = i4_min ( jlo + 4, m );
      cout << "\n";
      cout << "      ";
      for ( j = jlo; j <= jhi; j++ )
      {
        cout << setw(7) << j << "       ";
      }
      cout << "\n";
      cout << "\n";
      for ( i = 1; i <= l; i++ )
      {
        cout << setw(5) << i << ":";
        for ( j = jlo; j <= jhi; j++ )
        {
          cout << "  " << setw(12) << a[i-1+(j-1)*l+(k-1)*l*m];
        }
        cout << "\n";
      }
    }
  }

  return;
}


double *r8block_zeros_new ( int l, int m, int n )






























{
  double *a;
  int i;
  int j;
  int k;

  a = new double[l*m*n];

  for ( k = 0; k < n; k++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( i = 0; i < l; i++ )
      {
        a[i+j*l+k*l*m] = 0.0;
      }
    }
  }
  return a;
}


void r8cmat_delete ( int m, int n, double **a )


































{
  int j;

  for ( j = 0; j < n; j++ )
  {
    delete [] a[j];
  }

  delete [] a;

  return;
}


double **r8cmat_new ( int m, int n )
































{
  double **a;
  int j;

  a = new double *[n];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "R8CMAT_NEW - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    a[j] = new double[m];
    if ( a[j] == NULL )
    {
      cerr << "\n";
      cerr << "R8CMAT_NEW - Fatal error!\n";
      cerr << "  Unable to allocate row array.\n";
      exit ( 1 );
    }
  }

  return a;
}


void r8cmat_print ( int m, int n, double **a, string title )




































{
  r8cmat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}


void r8cmat_print_some ( int m, int n, double **a, int ilo, int jlo, int ihi,
  int jhi, string title )









































{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }



  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }
    cout << "\n";





    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j - 1 << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";



    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( ihi < m )
    {
      i2hi = ihi;
    }
    else
    {
      i2hi = m;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {



      cout << setw(5) << i - 1 << ": ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[j-1][i-1] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}


double *r8cmat_to_r8mat_new ( int m, int n, double **a )






































{
  double *b;
  int i;
  int j;

  b = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      b[i+j*m] = a[j][i];
    }
  }

  return b;
}


double **r8cmat_zeros_new ( int m, int n )
































{
  double **a;
  int i;
  int j;

  a = new double *[n];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "R8CMAT_ZEROS_NEW - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    a[j] = new double[m];
    if ( a[j] == NULL )
    {
      cerr << "\n";
      cerr << "R8CMAT_ZEROS_NEW - Fatal error!\n";
      cerr << "  Unable to allocate row array.\n";
      exit ( 1 );
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[j][i] = 0.0;
    }
  }
  return a;
}


double r8int_to_r8int ( double rmin, double rmax, double r, double r2min,
  double r2max )




































{
  double  r2;

  if ( rmax == rmin )
  {
    r2 = ( r2max + r2min ) / 2.0;
  }
  else
  {
    r2 = ( ( ( rmax - r        ) * r2min
           + (        r - rmin ) * r2max )
           / ( rmax     - rmin ) );
  }

  return r2;
}


int r8int_to_i4int ( double rmin, double rmax, double r, int imin, int imax )



































{
  int i;

  if ( rmax == rmin )
  {
    i = ( imax + imin ) / 2;
  }
  else
  {
    i = r8_nint (
      ( ( rmax - r        ) * ( double ) ( imin )
      + (        r - rmin ) * ( double ) ( imax ) )
      / ( rmax     - rmin ) );
  }

  return i;
}


void r8mat_add ( int m, int n, double alpha, double a[], double beta, 
  double b[], double c[] )





































{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+j*m] = alpha * a[i+j*m] + beta * b[i+j*m];
    }
  }
  return;
}


double *r8mat_add_new ( int m, int n, double alpha, double a[], double beta, 
  double b[] )





































{
  double *c;
  int i;
  int j;

  c = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+j*m] = alpha * a[i+j*m] + beta * b[i+j*m];
    }
  }
  return c;
}


double r8mat_amax ( int m, int n, double a[] )


































{
  int i;
  int j;
  double value;

  value = fabs ( a[0+0*m] );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = r8_max ( value, fabs ( a[i+j*m] ) );
    }
  }
  return value;
}


double *r8mat_border_add ( int m, int n, double table[] )















































{
  int i;
  int j;
  double *table2;

  table2 = new double[(m+2)*(n+2)];

  for ( j = 0; j < n+2; j++ )
  {
    for ( i = 0; i < m+2; i++ )
    {
      if ( i == 0 || i == m+1 || j == 0 || j == n+1 )
      {
        table2[i+j*(m+2)] = 0.0;
      }
      else
      {
        table2[i+j*(m+2)] = table[(i-1)+(j-1)*m];
      }
    }
  }
  return table2;
}


double *r8mat_border_cut ( int m, int n, double table[] )












































{
  int i;
  int j;
  double *table2;

  if ( m <= 2 || n <= 2 )
  {
    return NULL;
  }

  table2 = new double[(m-2)*(n-2)];

  for ( j = 0; j < n-2; j++ )
  {
    for ( i = 0; i < m-2; i++ )
    {
      table2[i+j*(m-2)] = table[(i+1)+(j+1)*m];
    }
  }

  return table2;
}


double *r8mat_cholesky_factor ( int n, double a[], int &flag )













































{
  double *c;
  int i;
  int j;
  int k;
  double sum2;
  double tol;

  flag = 0;
  tol = sqrt ( r8_epsilon ( ) );

  c = r8mat_copy_new ( n, n, a );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      c[i+j*n] = 0.0;
    }
    for ( i = j; i < n; i++ )
    {
      sum2 = c[j+i*n];
      for ( k = 0; k < j; k++ )
      {
        sum2 = sum2 - c[j+k*n] * c[i+k*n];
      }
      if ( i == j )
      {
        if ( 0.0 < sum2 )
        {
          c[i+j*n] = sqrt ( sum2 );
        }
        else if ( sum2 < - tol )
        {
          flag = 2;
          cerr << "\n";
          cerr << "R8MAT_CHOLESKY_FACTOR - Fatal error!\n";
          cerr << "  Matrix is not nonnegative definite.\n";
          cerr << "  Diagonal I = " << i << "\n";
          cerr << "  SUM2 = " << sum2 << "\n";
          exit ( 1 );
        }
        else
        {
          flag = 1;
          c[i+j*n] = 0.0;
        }
      }
      else
      {

        if ( c[j+j*n] != 0.0 )
        {
          c[i+j*n] = sum2 / c[j+j*n];
        }
        else
        {
          c[i+j*n] = 0.0;
        }
      }
    }
  }

  return c;
}


double *r8mat_cholesky_factor_upper ( int n, double a[], int &flag )

















































{
  double *c;
  int i;
  int j;
  int k;
  double sum2;

  flag = 0;

  c = r8mat_copy_new ( n, n, a );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      c[j+i*n] = 0.0;
    }
    for ( i = j; i < n; i++ )
    {
      sum2 = c[i+j*n];
      for ( k = 0; k < j; k++ )
      {
        sum2 = sum2 - c[k+j*n] * c[k+i*n];
      }
      if ( i == j )
      {
        if ( sum2 <= 0.0 )
        {
          flag = 1;
          return NULL;
        }
        c[j+i*n] = sqrt ( sum2 );
      }
      else
      {
        if ( c[j+j*n] != 0.0 )
        {
          c[j+i*n] = sum2 / c[j+j*n];
        }
        else
        {
          c[j+i*n] = 0.0;
        }
      }
    }
  }

  return c;
}


void r8mat_cholesky_inverse ( int n, double a[] )









































{
  int i;
  int j;
  int k;
  double s;
  double t;

  for ( j = 0; j < n; j++ )
  {
    s = 0.0;

    for ( k = 0; k < j; k++ )
    {
      t = a[k+j*n];
      for ( i = 0; i < k; i++ )
      {
        t = t - a[i+k*n] * a[i+j*n];
      }
      t = t / a[k+k*n];
      a[k+j*n] = t;
      s = s + t * t;
    }

    s = a[j+j*n] - s;

    if ( s <= 0.0 )
    {
      cerr << "\n";
      cerr << "R8MAT_CHOLESKY_INVERSE - Fatal error!\n";
      cerr << "  The matrix is singular.\n";
      exit ( 1 );
    }

    a[j+j*n] = sqrt ( s );

    for ( i = j + 1; i < n; i++ )
    {
      a[i+j*n] = 0.0;
    }
  }



  for ( k = 0; k < n; k++ )
  {
    a[k+k*n] = 1.0 / a[k+k*n];
    for ( i = 0; i < k; i++ )
    {
      a[i+k*n] = - a[i+k*n] * a[k+k*n];
    }

    for ( j = k + 1; j < n; j++ )
    {
      t = a[k+j*n];
      a[k+j*n] = 0.0;
      for ( i = 0; i <= k; i++ )
      {
        a[i+j*n] = a[i+j*n] + t * a[i+k*n];
      }
    }
  }



  for ( j = 0; j < n; j++ )
  {
    for ( k = 0; k < j; k++ )
    {
      t = a[k+j*n];
      for ( i = 0; i <= k; i++ )
      {
        a[i+k*n] = a[i+k*n] + t * a[i+j*n];
      }
    }
    t = a[j+j*n];
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*n] = a[i+j*n] * t;
    }
  }



  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      a[i+j*n] = a[j+i*n];
    }
  }

  return;
}


double *r8mat_cholesky_solve ( int n, double l[], double b[] )



































{
  double *x;
  double *y;



  y = r8mat_l_solve ( n, l, b );



  x = r8mat_lt_solve ( n, l, y );

  delete [] y;

  return x;
}


double *r8mat_cholesky_solve_upper ( int n, double r[], double b[] )



































{
  double *x;
  double *y;



  y = r8mat_ut_solve ( n, r, b );



  x = r8mat_u_solve ( n, r, y );

  delete [] y;

  return x;
}


void r8mat_copy ( int m, int n, double a1[], double a2[] )
































{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return;
}


double *r8mat_copy_new ( int m, int n, double a1[] )
































{
  double *a2;
  int i;
  int j;

  a2 = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return a2;
}


double *r8mat_covariance ( int m, int n, double x[] )



































{
  double *c;
  int i;
  int j;
  int k;
  double *x_mean;

  c = new double[m*m];
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+j*m] = 0.0;
    }
  }



  if ( n == 1 )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+i*m] = 1.0;
    }
    return c;
  }



  x_mean = new double[m];
  for ( i = 0; i < m; i++ )
  {
    x_mean[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      x_mean[i] = x_mean[i] + x[i+j*m];
    }
    x_mean[i] = x_mean[i] / ( double ) ( n );
  }



  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      for ( k = 0; k < n; k++ )
      {
        c[i+j*m] = c[i+j*m] 
          + ( x[i+k*m] - x_mean[i] ) * ( x[j+k*m] - x_mean[j] );
      }
    }
  }

  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+j*m] = c[i+j*m] / ( double ) ( n - 1 );
    }
  }

  delete [] x_mean;

  return c;
}


double r8mat_det ( int n, double a[] )








































{
  double *b;
  double det;
  int i;
  int j;
  int k;
  int kk;
  int m;
  double temp;

  b = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = a[i+j*n];
    }
  }

  det = 1.0;

  for ( k = 1; k <= n; k++ )
  {
    m = k;
    for ( kk = k+1; kk <= n; kk++ )
    {
      if ( fabs ( b[m-1+(k-1)*n] ) < fabs ( b[kk-1+(k-1)*n] ) )
      {
        m = kk;
      }
    }

    if ( m != k )
    {
      det = -det;

      temp = b[m-1+(k-1)*n];
      b[m-1+(k-1)*n] = b[k-1+(k-1)*n];
      b[k-1+(k-1)*n] = temp;
    }

    det = det * b[k-1+(k-1)*n];

    if ( b[k-1+(k-1)*n] != 0.0 )
    {
      for ( i = k+1; i <= n; i++ )
      {
        b[i-1+(k-1)*n] = -b[i-1+(k-1)*n] / b[k-1+(k-1)*n];
      }

      for ( j = k+1; j <= n; j++ )
      {
        if ( m != k )
        {
          temp = b[m-1+(j-1)*n];
          b[m-1+(j-1)*n] = b[k-1+(j-1)*n];
          b[k-1+(j-1)*n] = temp;
        }
        for ( i = k+1; i <= n; i++ )
        {
          b[i-1+(j-1)*n] = b[i-1+(j-1)*n] + b[i-1+(k-1)*n] * b[k-1+(j-1)*n];
        }
      }
    }
  }

  delete [] b;

  return det;
}


double r8mat_det_2d ( double a[] )




































{
  double det;

  det = a[0+0*2] * a[1+1*2] - a[0+1*2] * a[1+0*2];

  return det;
}


double r8mat_det_3d ( double a[] )




































{
  double det;

  det =
      a[0+0*3] * ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] )
    + a[0+1*3] * ( a[1+2*3] * a[2+0*3] - a[1+0*3] * a[2+2*3] )
    + a[0+2*3] * ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] );

  return det;
}


double r8mat_det_4d ( double a[] )






























{
  double det;

  det =
      a[0+0*4] * (
          a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] ) )
    - a[0+1*4] * (
          a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] ) )
    + a[0+2*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) )
    - a[0+3*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
        + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) );

  return det;
}


double r8mat_det_5d ( double a[] )






























{
  double b[4*4];
  double det;
  int i;
  int inc;
  int j;
  int k;
  double sign;




  det = 0.0;
  sign = 1.0;

  for ( k = 0; k < 5; k++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      for ( j = 0; j < 4; j++ )
      {
        if ( j < k )
        {
          inc = 0;
        }
        else
        {
          inc = 1;
        }
        b[i+j*4] = a[i+1+(j+inc)*5];
      }
    }

    det = det + sign * a[0+k*5] * r8mat_det_4d ( b );

    sign = - sign;
  }

  return det;
}


void r8mat_diag_add_scalar ( int n, double a[], double s )

































{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i+i*n] = a[i+i*n] + s;
  }

  return;
}


void r8mat_diag_add_vector ( int n, double a[], double v[] )
































{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i+i*n] = a[i+i*n] + v[i];
  }

  return;
}


void r8mat_diag_get_vector ( int n, double a[], double v[] )

































{
  int i;

  for ( i = 0; i < n; i++ )
  {
    v[i] = a[i+i*n];
  }

  return;
}


double *r8mat_diag_get_vector_new ( int n, double a[] )

































{
  int i;
  double *v;

  v = new double[n];

  for ( i = 0; i < n; i++ )
  {
    v[i] = a[i+i*n];
  }

  return v;
}


void r8mat_diag_set_scalar ( int n, double a[], double s )

































{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i+i*n] = s;
  }

  return;
}


void r8mat_diag_set_vector ( int n, double a[], double v[] )

































{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i+i*n] = v[i];
  }

  return;
}


double *r8mat_diagonal_new ( int n, double diag[] )
































{
  double *a;
  int i;
  int j;

  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = diag[i];
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }

  return a;
}


double r8mat_diff_frobenius ( int m, int n, double a[], double b[] )













































{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m] - b[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
}


double *r8mat_expand_linear ( int m, int n, double x[], int mfat, int nfat )













































{
  int i;
  int ihi;
  int ii;
  int iii;
  int ip1;
  int j;
  int jhi;
  int jj;
  int jjj;
  int jp1;
  int m2;
  int n2;
  double s;
  double t;
  double x00;
  double x01;
  double x10;
  double x11;
  double *xfat;

  m2 = ( m - 1 ) * ( mfat + 1 ) + 1;
  n2 = ( n - 1 ) * ( nfat + 1 ) + 1;

  xfat = new double[m2*n2];

  for ( i = 1; i <= m; i++ )
  {
    if ( i < m )
    {
      ihi = mfat;
    }
    else
    {
      ihi = 0;
    }

    for ( j = 1; j <= n; j++ )
    {
      if ( j < n )
      {
        jhi = nfat;
      }
      else
      {
        jhi = 0;
      }

      if ( i < m )
      {
        ip1 = i + 1;
      }
      else
      {
        ip1 = i;
      }

      if ( j < n )
      {
        jp1 = j + 1;
      }
      else
      {
        jp1 = j;
      }

      x00 = x[i-1+(j-1)*m];
      x10 = x[ip1-1+(j-1)*m];
      x01 = x[i-1+(jp1-1)*m];
      x11 = x[ip1-1+(jp1-1)*m];

      for ( ii = 0; ii <= ihi; ii++ )
      {
        s = ( double ) ( ii ) / ( double ) ( ihi + 1 );

        for ( jj = 0; jj <= jhi; jj++ )
        {
          t = ( double ) ( jj ) / ( double ) ( jhi + 1 );

          iii = 1 + ( i - 1 ) * ( mfat + 1 ) + ii;
          jjj = 1 + ( j - 1 ) * ( nfat + 1 ) + jj;

          xfat[iii-1+(jjj-1)*m2] =
                                            x00
              + s     * (       x10       - x00 )
              + t     * (             x01 - x00 )
              + s * t * ( x11 - x10 - x01 + x00 );
        }
      }
    }
  }

  return xfat;
}


double *r8mat_expand_linear2 ( int m, int n, double a[], int m2, int n2 )






































{
  double *a2;
  int i;
  int i1;
  int i2;
  int j;
  int j1;
  int j2;
  double r;
  double r1;
  double r2;
  double s;
  double s1;
  double s2;

  a2 = new double[m2*n2];

  for ( i = 1; i <= m2; i++ )
  {
    if ( m2 == 1 )
    {
      r = 0.5;
    }
    else
    {
      r = ( double ) ( i - 1 ) / ( double ) ( m2 - 1 );
    }

    i1 = 1 + ( int ) ( r * ( double ) ( m - 1 ) );
    i2 = i1 + 1;

    if ( m < i2 )
    {
      i1 = m - 1;
      i2 = m;
    }

    r1 = ( double ) ( i1 - 1 ) / ( double ) ( m - 1 );
    r2 = ( double ) ( i2 - 1 ) / ( double ) ( m - 1 );

    for ( j = 1; j <= n2; j++ )
    {
      if ( n2 == 1 )
      {
        s = 0.5;
      }
      else
      {
        s = ( double ) ( j - 1 ) / ( double ) ( n2 - 1 );
      }

      j1 = 1 + ( int ) ( s * ( double ) ( n - 1 ) );
      j2 = j1 + 1;

      if ( n < j2 )
      {
        j1 = n - 1;
        j2 = n;
      }

      s1 = ( double ) ( j1 - 1 ) / ( double ) ( n - 1 );
      s2 = ( double ) ( j2 - 1 ) / ( double ) ( n - 1 );

      a2[i-1+(j-1)*m2] =
        ( ( r2 - r ) * ( s2 - s ) * a[i1-1+(j1-1)*m]
        + ( r - r1 ) * ( s2 - s ) * a[i2-1+(j1-1)*m]
        + ( r2 - r ) * ( s - s1 ) * a[i1-1+(j2-1)*m]
        + ( r - r1 ) * ( s - s1 ) * a[i2-1+(j2-1)*m] )
        / ( ( r2 - r1 ) * ( s2 - s1 ) );
    }
  }

  return a2;
}


double *r8mat_flip_cols_new ( int m, int n, double a[] )
































{
  double *b;
  int i;
  int j;

  b = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      b[i+(n-1-j)*m] = a[i+j*m];
    }
  }

  return b;
}


double *r8mat_flip_rows_new ( int m, int n, double a[] )
































{
  double *b;
  int i;
  int j;

  b = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      b[(m-1-i)+j*m] = a[i+j*m];
    }
  }

  return b;
}


void r8mat_fs ( int n, double a[], double x[] )





































{
  double *a2;
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;

  a2 = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a2[i+j*n] = a[i+j*n];
    }
  }

  for ( jcol = 1; jcol <= n; jcol++ )
  {



    piv = fabs ( a2[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < fabs ( a2[i-1+(jcol-1)*n] ) )
      {
        piv = fabs ( a2[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cerr << "\n";
      cerr << "R8MAT_FS - Fatal error!\n";
      cerr << "  Zero pivot on step " << jcol << "\n";
      exit ( 1 );
    }



    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                  = a2[jcol-1+(j-1)*n];
        a2[jcol-1+(j-1)*n] = a2[ipiv-1+(j-1)*n];
        a2[ipiv-1+(j-1)*n] = t;
      }
      t         = x[jcol-1];
      x[jcol-1] = x[ipiv-1];
      x[ipiv-1] = t;
    }



    t = a2[jcol-1+(jcol-1)*n];
    a2[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a2[jcol-1+(j-1)*n] = a2[jcol-1+(j-1)*n] / t;
    }
    x[jcol-1] = x[jcol-1] / t;



    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a2[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a2[i-1+(jcol-1)*n];
        a2[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a2[i-1+(j-1)*n] = a2[i-1+(j-1)*n] + t * a2[jcol-1+(j-1)*n];
        }
        x[i-1] = x[i-1] + t * x[jcol-1];
      }
    }
  }



  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      x[i-1] = x[i-1] - a2[i-1+(jcol-1)*n] * x[jcol-1];
    }
  }

  delete [] a2;

  return;
}


double *r8mat_fs_new ( int n, double a[], double b[] )









































{
  double *a2;
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;
  double *x;

  a2 = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a2[i+j*n] = a[i+j*n];
    }
  }

  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  for ( jcol = 1; jcol <= n; jcol++ )
  {



    piv = fabs ( a2[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < fabs ( a2[i-1+(jcol-1)*n] ) )
      {
        piv = fabs ( a2[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cerr << "\n";
      cerr << "R8MAT_FS_NEW - Fatal error!\n";
      cerr << "  Zero pivot on step " << jcol << "\n";
      exit ( 1 );
    }



    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                  = a2[jcol-1+(j-1)*n];
        a2[jcol-1+(j-1)*n] = a2[ipiv-1+(j-1)*n];
        a2[ipiv-1+(j-1)*n] = t;
      }
      t         = x[jcol-1];
      x[jcol-1] = x[ipiv-1];
      x[ipiv-1] = t;
    }



    t = a2[jcol-1+(jcol-1)*n];
    a2[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a2[jcol-1+(j-1)*n] = a2[jcol-1+(j-1)*n] / t;
    }
    x[jcol-1] = x[jcol-1] / t;



    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a2[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a2[i-1+(jcol-1)*n];
        a2[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a2[i-1+(j-1)*n] = a2[i-1+(j-1)*n] + t * a2[jcol-1+(j-1)*n];
        }
        x[i-1] = x[i-1] + t * x[jcol-1];
      }
    }
  }



  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      x[i-1] = x[i-1] - a2[i-1+(jcol-1)*n] * x[jcol-1];
    }
  }

  delete [] a2;

  return x;
}


void r8mat_fss ( int n, double a[], int nb, double x[] )







































{
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;

  for ( jcol = 1; jcol <= n; jcol++ )
  {



    piv = fabs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < fabs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = fabs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cerr << "\n";
      cerr << "R8MAT_FSS - Fatal error!\n";
      cerr << "  Zero pivot on step " << jcol << "\n";
      exit ( 1 );
    }



    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      for ( j = 0; j < nb; j++ )
      {
        t            = x[jcol-1+j*n];
        x[jcol-1+j*n] = x[ipiv-1+j*n];
        x[ipiv-1+j*n] = t;
      }
    }



    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    for ( j = 0; j < nb; j++ )
    {
      x[jcol-1+j*n] = x[jcol-1+j*n] / t;
    }



    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        for ( j = 0; j < nb; j++ )
        {
          x[i-1+j*n] = x[i-1+j*n] + t * x[jcol-1+j*n];
        }
      }
    }
  }



  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      for ( j = 0; j < nb; j++ )
      {
        x[i-1+j*n] = x[i-1+j*n] - a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
      }
    }
  }

  return;
}


double *r8mat_fss_new ( int n, double a[], int nb, double b[] )








































{
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;
  double *x;

  x = new double[n*nb];

  for ( j = 0; j < nb; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = b[i+j*n];
    }
  }
  for ( jcol = 1; jcol <= n; jcol++ )
  {



    piv = fabs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol + 1; i <= n; i++ )
    {
      if ( piv < fabs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = fabs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cerr << "\n";
      cerr << "R8MAT_FSS_NEW - Fatal error!\n";
      cerr << "  Zero pivot on step " << jcol << "\n";
      exit ( 1 );
    }



    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      for ( j = 0; j < nb; j++ )
      {
        t            = x[jcol-1+j*n];
        x[jcol-1+j*n] = x[ipiv-1+j*n];
        x[ipiv-1+j*n] = t;
      }
    }



    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    for ( j = 0; j < nb; j++ )
    {
      x[jcol-1+j*n] = x[jcol-1+j*n] / t;
    }



    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        for ( j = 0; j < nb; j++ )
        {
          x[i-1+j*n] = x[i-1+j*n] + t * x[jcol-1+j*n];
        }
      }
    }
  }



  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      for ( j = 0; j < nb; j++ )
      {
        x[i-1+j*n] = x[i-1+j*n] - a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
      }
    }
  }

  return x;
}


double *r8mat_givens_post ( int n, double a[], int row, int col )








































{
  double *g;
  double theta;

  g = r8mat_identity_new ( n );

  theta = atan2 ( a[row-1+(col-1)*n], a[row-1+(row-1)*n] );

  g[row-1+(row-1)*n] =  cos ( theta );
  g[row-1+(col-1)*n] = -sin ( theta );
  g[col-1+(row-1)*n] =  sin ( theta );
  g[col-1+(col-1)*n] =  cos ( theta );

  return g;
}


double *r8mat_givens_pre ( int n, double a[], int row, int col )








































{
  double *g;
  double theta;

  g = r8mat_identity_new ( n );

  theta = atan2 ( a[row-1+(col-1)*n], a[col-1+(col-1)*n] );

  g[row-1+(row-1)*n] =  cos ( theta );
  g[row-1+(col-1)*n] = -sin ( theta );
  g[col-1+(row-1)*n] =  sin ( theta );
  g[col-1+(col-1)*n] =  cos ( theta );

  return g;
}


double *r8mat_hess ( double (*fx) ( int n, double x[] ), int n, double x[] )












































{
  double eps;
  double f00;
  double fmm;
  double fmp;
  double fpm;
  double fpp;
  double *h;
  int i;
  int j;
  double *s;
  double xi;
  double xj;



  s = new double[n];

  eps = pow ( r8_epsilon ( ), 0.33 );

  for ( i = 0; i < n; i++ )
  {
    s[i] = eps * r8_max ( fabs ( x[i] ), 1.0 );
  }



  h = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    xi = x[i];

    f00 = fx ( n, x );

    x[i] = xi + s[i];
    fpp = fx ( n, x );

    x[i] = xi - s[i];
    fmm = fx ( n, x );

    h[i+i*n] = ( ( fpp - f00 ) + ( fmm - f00 ) ) / s[i] / s[i];

    x[i] = xi;
  }



  for ( i = 0; i < n; i++ )
  {
    xi = x[i];

    for ( j = i+1; j < n; j++ )
    {
      xj = x[j];

      x[i] = xi + s[i];
      x[j] = xj + s[j];
      fpp = fx ( n, x );

      x[i] = xi + s[i];
      x[j] = xj - s[j];
      fpm = fx ( n, x );

      x[i] = xi - s[i];
      x[j] = xj + s[j];
      fmp = fx ( n, x );

      x[i] = xi - s[i];
      x[j] = xj - s[j];
      fmm = fx ( n, x );

      h[j+i*n] = ( ( fpp - fpm ) + ( fmm - fmp ) ) / ( 4.0 * s[i] * s[j] );

      h[i+j*n] = h[j+i*n];

      x[j] = xj;
    }
    x[i] = xi;
  }

  delete [] s;

  return h;
}


void r8mat_house_axh ( int n, double a[], double v[] )







































{
  double *ah;
  int i;
  int j;
  int k;
  double v_normsq;

  v_normsq = 0.0;
  for ( i = 0; i < n; i++ )
  {
    v_normsq = v_normsq + v[i] * v[i];
  }



  ah = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      ah[i+j*n] = a[i+j*n];
      for ( k = 0; k < n; k++ )
      {
        ah[i+j*n] = ah[i+j*n] - 2.0 * a[i+k*n] * v[k] * v[j] / v_normsq;
      }
    }
  }



  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = ah[i+j*n];
    }
  }
  delete [] ah;

  return;
}


double *r8mat_house_axh_new ( int n, double a[], double v[] )








































{
  double *ah;
  int i;
  int j;
  int k;
  double v_normsq;

  v_normsq = 0.0;
  for ( i = 0; i < n; i++ )
  {
    v_normsq = v_normsq + v[i] * v[i];
  }



  ah = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      ah[i+j*n] = a[i+j*n];
      for ( k = 0; k < n; k++ )
      {
        ah[i+j*n] = ah[i+j*n] - 2.0 * a[i+k*n] * v[k] * v[j] / v_normsq;
      }
    }
  }

  return ah;
}


double *r8mat_house_form ( int n, double v[] )


































{
  double beta;
  double *h;
  int i;
  int j;



  beta = 0.0;
  for ( i = 0; i < n; i++ )
  {
    beta = beta + v[i] * v[i];
  }



  h = r8mat_identity_new ( n );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      h[i+j*n] = h[i+j*n] - 2.0 * v[i] * v[j] / beta;
    }
  }

  return h;
}


double *r8mat_house_hxa ( int n, double a[], double v[] )








































{
  double *ha;
  int i;
  int j;
  int k;
  double v_normsq;

  v_normsq = 0.0;
  for ( i = 0; i < n; i++ )
  {
    v_normsq = v_normsq + v[i] * v[i];
  }



  ha = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      ha[i+j*n] = a[i+j*n];
      for ( k = 0; k < n; k++ )
      {
        ha[i+j*n] = ha[i+j*n] - 2.0 * v[i] * v[k] * a[k+j*n] / v_normsq;
      }
    }
  }

  return ha;
}


double *r8mat_house_post ( int n, double a[], int row, int col )












































{
  double *a_row;
  double *h;
  int j;
  double *v;



  a_row = new double[n];

  for ( j = 0; j < col-1; j++ )
  {
    a_row[j] = 0.0;
  }
  for ( j = col - 1; j < n; j++ )
  {
    a_row[j] = a[row+j*n];
  }



  v = r8vec_house_column ( n, a_row, col );



  h = r8mat_house_form ( n, v );



  delete [] a_row;
  delete [] v;

  return h;
}


double *r8mat_house_pre ( int n, double a[], int row, int col )












































{
  double *a_col;
  double *h;
  int i;
  double *v;



  a_col = new double[n];

  for ( i = 0; i < row-1; i++ )
  {
    a_col[i] = 0.0;
  }
  for ( i = row-1; i < n; i++ )
  {
    a_col[i] = a[i+col*n];
  }



  v = r8vec_house_column ( n, a_col, row );



  h = r8mat_house_form ( n, v );



  delete [] a_col;
  delete [] v;

  return h;
}


void r8mat_identity ( int n, double a[] )






























{
  int i;
  int j;
  int k;

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[k] = 1.0;
      }
      else
      {
        a[k] = 0.0;
      }
      k = k + 1;
    }
  }

  return;
}


double *r8mat_identity_new ( int n )






























{
  double *a;
  int i;
  int j;
  int k;

  a = new double[n*n];

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[k] = 1.0;
      }
      else
      {
        a[k] = 0.0;
      }
      k = k + 1;
    }
  }

  return a;
}


double *r8mat_indicator_new ( int m, int n )








































{
  double *a;
  int fac;
  int i;
  int j;

  a = new double[m*n];

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      a[i-1+(j-1)*m] = ( double ) ( fac * i + j );
    }
  }
  return a;
}


double *r8mat_inverse_2d ( double a[] )






























{
  double *b;
  double det;



  det = a[0+0*2] * a[1+1*2] - a[0+1*2] * a[1+0*2];



  if ( det == 0.0 )
  {
    return NULL;
  }



  b = new double[2*2];

  b[0+0*2] = + a[1+1*2] / det;
  b[0+1*2] = - a[0+1*2] / det;
  b[1+0*2] = - a[1+0*2] / det;
  b[1+1*2] = + a[0+0*2] / det;

  return b;
}


double *r8mat_inverse_3d ( double a[] )




































{
  double *b;
  double det;



  det =
     a[0+0*3] * ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] )
   + a[0+1*3] * ( a[1+2*3] * a[2+0*3] - a[1+0*3] * a[2+2*3] )
   + a[0+2*3] * ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] );

  if ( det == 0.0 )
  {
    return NULL;
  }

  b = new double[3*3];

  b[0+0*3] =   ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] ) / det;
  b[0+1*3] = - ( a[0+1*3] * a[2+2*3] - a[0+2*3] * a[2+1*3] ) / det;
  b[0+2*3] =   ( a[0+1*3] * a[1+2*3] - a[0+2*3] * a[1+1*3] ) / det;

  b[1+0*3] = - ( a[1+0*3] * a[2+2*3] - a[1+2*3] * a[2+0*3] ) / det;
  b[1+1*3] =   ( a[0+0*3] * a[2+2*3] - a[0+2*3] * a[2+0*3] ) / det;
  b[1+2*3] = - ( a[0+0*3] * a[1+2*3] - a[0+2*3] * a[1+0*3] ) / det;

  b[2+0*3] =   ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] ) / det;
  b[2+1*3] = - ( a[0+0*3] * a[2+1*3] - a[0+1*3] * a[2+0*3] ) / det;
  b[2+2*3] =   ( a[0+0*3] * a[1+1*3] - a[0+1*3] * a[1+0*3] ) / det;

  return b;
}


double *r8mat_inverse_4d ( double a[] )






























{
  double *b;
  double det;



  det = r8mat_det_4d ( a );



  if ( det == 0.0 )
  {
    return NULL;
  }



  b = new double[4*4];

  b[0+0*4] =
    +(
    + a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
    + a[1+2*4] * ( a[2+3*4] * a[3+1*4] - a[2+1*4] * a[3+3*4] )
    + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
    ) / det;

  b[1+0*4] =
    -(
    + a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
    + a[1+2*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
    + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
    ) / det;

  b[2+0*4] =
    +(
    + a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
    + a[1+1*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
    + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
    ) / det;

  b[3+0*4] =
    -(
    + a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
    + a[1+1*4] * ( a[2+2*4] * a[3+0*4] - a[2+0*4] * a[3+2*4] )
    + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
    ) / det;

  b[0+1*4] =
    -(
    + a[0+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
    + a[0+2*4] * ( a[2+3*4] * a[3+1*4] - a[2+1*4] * a[3+3*4] )
    + a[0+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
    ) / det;

  b[1+1*4] =
    +(
    + a[0+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
    + a[0+2*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
    + a[0+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
    ) / det;

  b[2+1*4] =
    -(
    + a[0+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
    + a[0+1*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
    + a[0+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
    ) / det;

  b[3+1*4] =
    +(
    + a[0+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
    + a[0+1*4] * ( a[2+2*4] * a[3+0*4] - a[2+0*4] * a[3+2*4] )
    + a[0+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
    ) / det;

  b[0+2*4] =
    +(
    + a[0+1*4] * ( a[1+2*4] * a[3+3*4] - a[1+3*4] * a[3+2*4] )
    + a[0+2*4] * ( a[1+3*4] * a[3+1*4] - a[1+1*4] * a[3+3*4] )
    + a[0+3*4] * ( a[1+1*4] * a[3+2*4] - a[1+2*4] * a[3+1*4] )
    ) / det;

  b[1+2*4] =
    -(
    + a[0+0*4] * ( a[1+2*4] * a[3+3*4] - a[1+3*4] * a[3+2*4] )
    + a[0+2*4] * ( a[1+3*4] * a[3+0*4] - a[1+0*4] * a[3+3*4] )
    + a[0+3*4] * ( a[1+0*4] * a[3+2*4] - a[1+2*4] * a[3+0*4] )
    ) / det;

  b[2+2*4] =
    +(
    + a[0+0*4] * ( a[1+1*4] * a[3+3*4] - a[1+3*4] * a[3+1*4] )
    + a[0+1*4] * ( a[1+3*4] * a[3+0*4] - a[1+0*4] * a[3+3*4] )
    + a[0+3*4] * ( a[1+0*4] * a[3+1*4] - a[1+1*4] * a[3+0*4] )
    ) / det;

  b[3+2*4] =
    -(
    + a[0+0*4] * ( a[1+1*4] * a[3+2*4] - a[1+2*4] * a[3+1*4] )
    + a[0+1*4] * ( a[1+2*4] * a[3+0*4] - a[1+0*4] * a[3+2*4] )
    + a[0+2*4] * ( a[1+0*4] * a[3+1*4] - a[1+1*4] * a[3+0*4] )
    ) / det;

  b[0+3*4] =
    -(
    + a[0+1*4] * ( a[1+2*4] * a[2+3*4] - a[1+3*4] * a[2+2*4] )
    + a[0+2*4] * ( a[1+3*4] * a[2+1*4] - a[1+1*4] * a[2+3*4] )
    + a[0+3*4] * ( a[1+1*4] * a[2+2*4] - a[1+2*4] * a[2+1*4] )
    ) / det;

  b[1+3*4] =
    +(
    + a[0+0*4] * ( a[1+2*4] * a[2+3*4] - a[1+3*4] * a[2+2*4] )
    + a[0+2*4] * ( a[1+3*4] * a[2+0*4] - a[1+0*4] * a[2+3*4] )
    + a[0+3*4] * ( a[1+0*4] * a[2+2*4] - a[1+2*4] * a[2+0*4] )
    ) / det;

  b[2+3*4] =
    -(
    + a[0+0*4] * ( a[1+1*4] * a[2+3*4] - a[1+3*4] * a[2+1*4] )
    + a[0+1*4] * ( a[1+3*4] * a[2+0*4] - a[1+0*4] * a[2+3*4] )
    + a[0+3*4] * ( a[1+0*4] * a[2+1*4] - a[1+1*4] * a[2+0*4] )
    ) / det;

  b[3+3*4] =
    +(
    + a[0+0*4] * ( a[1+1*4] * a[2+2*4] - a[1+2*4] * a[2+1*4] )
    + a[0+1*4] * ( a[1+2*4] * a[2+0*4] - a[1+0*4] * a[2+2*4] )
    + a[0+2*4] * ( a[1+0*4] * a[2+1*4] - a[1+1*4] * a[2+0*4] )
    ) / det;

  return b;
}


bool r8mat_is_binary ( int m, int n, double x[] )
































{
  int i;
  int j;
  bool value;

  value = true;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( x[i+j*m] != 0.0 && x[i+j*m] != 1.0 )
      {
        value = false;
        break;
      }
    }
  }
  return value;
}


double r8mat_is_identity ( int n, double a[] )



































{
  double error_frobenius;
  int i;
  int j;
  double t;

  error_frobenius = 0.0;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i == j )
      {
        t = a[i+j*n] - 1.0;
      }
      else
      {
        t = a[i+j*n];
      }
      error_frobenius = error_frobenius + t * t;
    }
  }
  error_frobenius = sqrt ( error_frobenius );

  return error_frobenius;
}


bool r8mat_is_in_01 ( int m, int n, double a[] )



































{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] < 0.0 || 1.0 < a[i+j*m] )
      {
        return false;
      }
    }
  }

  return true;
}


bool r8mat_is_insignificant ( int m, int n, double r[], double s[] )






























{
  int i;
  int j;
  double t;
  double tol;
  bool value;

  value = true;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      t = r[i+j*m] + s[i+j*m];
      tol = r8_epsilon ( ) * fabs ( r[i+j*m] );

      if ( tol < fabs ( r[i+j*m] - t ) )
      {
        value = false;
        break;
      }
    }
  }
  return value;
}


bool r8mat_is_significant ( int m, int n, double r[], double s[] )






























{
  int i;
  int j;
  double t;
  double tol;
  bool value;

  value = false;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      t = r[i+j*m] + s[i+j*m];
      tol = r8_epsilon ( ) * fabs ( r[i+j*m] );

      if ( tol < fabs ( r[i+j*m] - t ) )
      {
        value = true;
        break;
      }
    }
  }
  return value;
}


double r8mat_is_symmetric ( int m, int n, double a[] )

































{
  int i;
  int j;
  const double r8_huge = 1.79769313486231571E+308;
  double value;

  if ( m != n )
  {
    value = r8_huge;
    return value;
  }

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m] - a[j+i*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
}


double *r8mat_jac ( int m, int n, double eps,
  double *(*fx) ( int m, int n, double x[] ), double x[] )




























































{
  double del;
  double *fprime;
  int i;
  int j;
  double xsave;
  double *work1;
  double *work2;

  fprime = new double[m*n];



  work2 = fx ( m, n, x );




  for ( j = 0; j < n; j++ )
  {
    xsave = x[j];
    del = eps * ( 1.0 + fabs ( x[j] ) );
    x[j] = x[j] + del;
    work1 = fx ( m, n, x );
    x[j] = xsave;
    for ( i = 0; i < m; i++ )
    {
      fprime[i+j*m] = ( work1[i] - work2[i] ) / del;
    }
    delete [] work1;
  }
  delete [] work2;

  return fprime;
}


double *r8mat_kronecker ( int m1, int n1, double a[], int m2, int n2, 
  double b[] )












































{
  double *c;
  int i;
  int i0;
  int i1;
  int i2;
  int j;
  int j0;
  int j1;
  int j2;
  int m;
  int n;

  m = m1 * m2;
  n = n1 * n2;
  c = new double[m*n];

  for ( j1 = 0; j1 < n1; j1++ )
  {
    for ( i1 = 0; i1 < m1; i1++ )
    {
      i0 = i1 * m2;
      j0 = j1 * n2;
      j = j0;
      for ( j2 = 0; j2 < n2; j2++ )
      {
        i = i0;
        for ( i2 = 0; i2 < m2; i2++ )
        {
          c[i+j*m] = a[i1+j1*m1] * b[i2+j2*m2];
          i = i + 1;
        }
        j = j + 1;
      }
    }
  }

  return c;
}


double *r8mat_l_inverse ( int n, double a[] )












































{
  double *b;
  int i;
  int j;
  int k;
  double temp;

  b = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i < j )
      {
        b[i+j*n] = 0.0;
      }
      else if ( j == i )
      {
        b[i+j*n] = 1.0 / a[i+j*n];
      }
      else
      {
        temp = 0.0;
        for ( k = 0; k < i; k++ )
        {
          temp = temp + a[i+k*n] * b[k+j*n];
        }
        b[i+j*n] = -temp / a[i+i*n];
      }
    }
  }

  return b;
}


void r8mat_l_print ( int m, int n, double a[], string title )














































{
  int i;
  int indx[10];
  int j;
  int jhi;
  int jlo;
  int jmax;
  int nn;
  int size;

  cout << "\n";
  cout << title << "\n";

  jmax = i4_min ( n, m );

  if ( m <= n )
  {
    size = ( m * ( m + 1 ) ) / 2;
  }
  else
  {
    size = ( n * ( n + 1 ) ) / 2 + ( m - n ) * n;
  }

  if ( r8vec_is_integer ( size, a ) )
  {
    nn = 10;
    for ( jlo = 1; jlo <= jmax; jlo = jlo + nn )
    {
      jhi = i4_min ( jlo + nn - 1, i4_min ( m, jmax ) );
      cout << "\n";
      cout << "  Col   ";
      for ( j = jlo; j <= jhi; j++ )
      {
        cout << setw(6) << j;
      }
      cout << "\n";
      cout << "  Row  \n";
      for ( i = jlo; i <= m; i++ )
      {
        jhi = i4_min ( jlo + nn - 1, i4_min ( i, jmax ) );
        for ( j = jlo; j <= jhi; j++ )
        {
          indx[j-jlo] = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2;
        }
        cout << "  " << setw(6) << i;
        for ( j = 0; j <= jhi-jlo; j++ )
        {
          cout << setw(6) << a[indx[j]-1];
        }
        cout << "\n";
      }
    }
  }
  else if ( r8vec_amax ( size, a ) < 1000000.0 )
  {
    nn = 5;
    for ( jlo = 1; jlo <= jmax; jlo = jlo + nn )
    {
      jhi = i4_min ( jlo + nn - 1, i4_min ( m - 1, jmax ) );
      cout << "\n";
      cout << "  Col ";
      for ( j = jlo; j <= jhi; j++ )
      {
        cout << setw(14) << j;
      }
      cout << "\n";
      cout << "  Row  \n";
      for ( i = jlo; i <= m; i++ )
      {
        jhi = i4_min ( jlo + nn - 1, i4_min ( i, jmax ) );
        for ( j = jlo; j <= jhi; j++ )
        {
          indx[j-jlo] = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2;
        }
        cout << "  " << setw(6) << i;
        for ( j = 0; j <= jhi-jlo; j++ )
        {
          cout << setw(14) << a[indx[j]-1];
        }
        cout << "\n";
      }
    }
  }
  else
  {
    nn = 5;

    for ( jlo = 1; jlo <= jmax; jlo = jlo + nn )
    {
      jhi = i4_min ( jlo + nn - 1, i4_min ( m - 1, jmax ) );
      cout << "\n";
      cout << "  Col ";
      for ( j = jlo; j <= jhi; j++ )
      {
        cout << setw(7) << j << "       ";
      }
      cout << "\n";
      cout << "  Row \n";
      for ( i = jlo; i <= m; i++ )
      {
        jhi = i4_min ( jlo + nn - 1, i4_min ( i, jmax ) );
        for ( j = jlo; j <= jhi; j++ )
        {
          indx[j-jlo] = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2;
        }
        cout << setw(6) << i;
        for ( j = 0; j <= jhi-jlo; j++ )
        {
          cout << setw(14) << a[indx[j]-1];
        }
      }
    }
  }

  return;
}


double *r8mat_l_solve ( int n, double a[], double b[] )


































{
  int i;
  int j;
  double temp;
  double *x;

  x = new double[n];



  for ( i = 0; i < n; i++ )
  {
    temp = 0.0;
    for ( j = 0; j < i; j++ )
    {
      temp = temp + a[i+j*n] * x[j];
    }
    x[i] = ( b[i] - temp ) / a[i+i*n];
  }

  return x;
}


double *r8mat_l1_inverse ( int n, double a[] )













































{
  double *b;
  int i;
  int j;
  int k;

  b = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i < j )
      {
        b[i+j*n] = 0.0;
      }
      else if ( j == i )
      {
        b[i+j*n] = 1.0;
      }
      else
      {
        b[i+j*n] = 0.0;
        for ( k = 0; k < i; k++ )
        {
          b[i+j*n] = b[i+j*n] - a[i+k*n] * b[k+j*n];
        }
      }
    }
  }

  return b;
}


double *r8mat_lt_solve ( int n, double a[], double b[] )






































{
  int i;
  int j;
  double *x;

  x = new double[n];

  for ( j = n-1; 0 <= j; j-- )
  {
    x[j] = b[j];
    for ( i = j+1; i < n; i++ )
    {
      x[j] = x[j] - x[i] * a[i+j*n];
    }
    x[j] = x[j] / a[j+j*n];
  }

  return x;
}


void r8mat_lu ( int m, int n, double a[], double l[], double p[], double u[] )
















































{
  int i;
  int ipiv;
  int j;
  int k;
  double pivot;
  double t;







  r8mat_copy ( m, n, a, u );

  r8mat_zeros ( m, m, l );
  r8mat_zeros ( m, m, p );
  for ( i = 0; i < m; i++ )
  {
    l[i+i*m] = 1.0;
    p[i+i*m] = 1.0;
  }



  for ( j = 0; j < i4_min ( m - 1, n ); j++ )
  {
    pivot = 0.0;
    ipiv = -1;

    for ( i = j; i < m; i++ )
    {
      if ( pivot < fabs ( u[i+j*m] ) )
      {
        pivot = fabs ( u[i+j*m] );
        ipiv = i;
      }
    }



    if ( ipiv != -1 )
    {
      for ( k = 0; k < n; k++ )
      {
        t = u[j+k*m];
        u[j+k*m] = u[ipiv+j*m];
        u[ipiv+k*m] = t;

        t = l[j+k*m];
        l[j+k*m] = l[ipiv+j*m];
        l[ipiv+k*m] = t;

        t = p[j+k*m];
        p[j+k*m] = p[ipiv+j*m];
        p[ipiv+k*m] = t;
      }



      for ( i = j+1; i < m; i++ )
      {
        if ( u[i+j*m] != 0.0 )
        {
          l[i+j*m] = u[i+j*m] / u[j+j*m];

          u[i+j*m] = 0.0;

          for ( k = j+1; k < n; k++ )
          {
            u[i+k*m] = u[i+k*m] - l[i+j*m] * u[j+k*m];
          }
        }
      }
    }
  }

  return;
}


double r8mat_max ( int m, int n, double a[] )


































{
  int i;
  int j;
  double value;

  value = a[0+0*m];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( value < a[i+j*m] )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}


void r8mat_max_index ( int m, int n, double a[], int &i_max, int &j_max )


































{
  int i;
  int i2;
  int j;
  int j2;

  i2 = -1;
  j2 = -1;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( i2 == -1 && j2 == -1 )
      {
        i2 = i;
        j2 = j;
      }
      else if ( a[i2+j2*m] < a[i+j*m] )
      {
        i2 = i;
        j2 = j;
      }
    }
  }

  i_max = i2 + 1;
  j_max = j2 + 1;

  return;
}


double r8mat_maxcol_minrow ( int m, int n, double a[] )







































{
  int i;
  int j;
  double minrow;
  const double r8_huge = 1.79769313486231571E+308;
  double value;

  value = - r8_huge;

  for ( i = 0; i < m; i++ )
  {
    minrow = r8_huge;

    for ( j = 0; j < n; j++ )
    {
      minrow = r8_min ( minrow, a[i+j*m] );
    }
    value = r8_max ( value, minrow );
  }

  return value;
}


double r8mat_maxrow_mincol ( int m, int n, double a[] )







































{
  int i;
  int j;
  double mincol;
  const double r8_huge = 1.79769313486231571E+308;
  double value;

  value = - r8_huge;

  for ( j = 0; j < n; j++ )
  {
    mincol = r8_huge;
    for ( i = 0; i < m; i++ )
    {
      mincol = r8_min ( mincol, a[i+j*m] );
    }
    value = r8_max ( value, mincol );
  }
  return value;
}


double r8mat_mean ( int m, int n, double a[] )


































{
  int i;
  int j;
  double value;

  value = 0.0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + a[i+j*m];
    }
  }
  value = value / ( double ) ( m * n );

  return value;
}


double r8mat_min ( int m, int n, double a[] )


































{
  int i;
  int j;
  double value;

  value = a[0+0*m];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] < value )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}


void r8mat_min_index ( int m, int n, double a[], int &i_min, int &j_min )


































{
  int i;
  int i2;
  int j;
  int j2;

  i2 = -1;
  j2 = -1;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( i2 == -1 && j2 == -1 )
      {
        i2 = i;
        j2 = j;
      }
      else if ( a[i+j*m] < a[i2+j2*m] )
      {
        i2 = i;
        j2 = j;
      }
    }
  }

  i_min = i2 + 1;
  j_min = j2 + 1;

  return;
}


double r8mat_mincol_maxrow ( int m, int n, double a[] )







































{
  int i;
  int j;
  double maxrow;
  const double r8_huge = 1.79769313486231571E+308;
  double value;

  value = r8_huge;

  for ( i = 0; i < m; i++ )
  {
    maxrow = - r8_huge;
    for ( j = 0; j < n; j++ )
    {
      maxrow = r8_max ( maxrow, a[i+j*m] );
    }
    value = r8_min ( value, maxrow );
  }

  return value;
}


double r8mat_minrow_maxcol ( int m, int n, double a[] )







































{
  int i;
  int j;
  double maxcol;
  const double r8_huge = 1.79769313486231571E+308;
  double value;

  value = r8_huge;

  for ( j = 0; j < n; j++ )
  {
    maxcol = - r8_huge;
    for ( i = 0; i < m; i++ )
    {
      maxcol = r8_max ( maxcol, a[i+j*m] );
    }
    value = r8_min ( value, maxcol );
  }

  return value;
}


void r8mat_minvm ( int n1, int n2, double a[], double b[], double c[] )































{
  double *alu;
  double *d;

  alu = r8mat_copy_new ( n1, n1, a );

  d = r8mat_fss_new ( n1, alu, n2, b );

  r8mat_copy ( n1, n2, d, c );

  delete [] alu;
  delete [] d;

  return;
}


double *r8mat_minvm_new ( int n1, int n2, double a[], double b[] )































{
  double *alu;
  double *c;

  alu = r8mat_copy_new ( n1, n1, a );
  c = r8mat_fss_new ( n1, alu, n2, b );
 
  delete [] alu;

  return c;
}


void r8mat_mm ( int n1, int n2, int n3, double a[], double b[], double c[] )


































{
  double *c1;
  int i;
  int j;
  int k;

  c1 = new double[n1*n3];

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c1[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c1[i+j*n1] = c1[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  r8mat_copy ( n1, n3, c1, c );

  delete [] c1;

  return;
}


double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] )


































{
  double *c;
  int i;
  int j;
  int k;

  c = new double[n1*n3];

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return c;
}


double *r8mat_mmt_new ( int n1, int n2, int n3, double a[], double b[] )


































{
  double *c;
  int i;
  int j;
  int k;

  c = new double[n1*n3];

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[j+k*n3];
      }
    }
  }

  return c;
}


double *r8mat_mtm_new ( int n1, int n2, int n3, double a[], double b[] )


































{
  double *c;
  int i;
  int j;
  int k;

  c = new double[n1*n3];

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[k+i*n2] * b[k+j*n2];
      }
    }
  }

  return c;
}


void r8mat_mtv ( int m, int n, double a[], double x[], double atx[] )




































{
  int i;
  int j;
  double *y;

  y = new double[n];

  for ( j = 0; j < n; j++ )
  {
    y[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      y[j] = y[j] + a[i+j*m] * x[i];
    }
  }

  r8vec_copy ( n, y, atx );

  free ( y );

  return;
}


double *r8mat_mtv_new ( int m, int n, double a[], double x[] )




































{
  int i;
  int j;
  double *y;

  y = new double[n];

  for ( j = 0; j < n; j++ )
  {
    y[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      y[j] = y[j] + a[i+j*m] * x[i];
    }
  }

  return y;
}


void r8mat_mv ( int m, int n, double a[], double x[], double ax[] )




































{
  int i;
  int j;
  double *y;

  y = ( double * ) malloc ( m * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  r8vec_copy ( m, y, ax );

  free ( y );

  return;
}


double *r8mat_mv_new ( int m, int n, double a[], double x[] )




































{
  int i;
  int j;
  double *y;

  y = new double[m];

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  return y;
}


void r8mat_nint ( int m, int n, double a[] )






























{
  int i;
  int j;
  int s;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] < 0.0 )
      {
        s = -1;
      }
      else
      {
        s = 1;
      }
      a[i+j*m] = s * ( int ) ( fabs ( a[i+j*m] ) + 0.5 );
    }
  }

  return;
}


int r8mat_nonzeros ( int m, int n, double a[] )
































{
  int i;
  int j;
  int value;

  value = 0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] != 0.0 )
      {
        value = value + 1;
      }
    }
  }

  return value;
}


double r8mat_norm_eis ( int m, int n, double a[] )







































{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + fabs ( a[i+j*m] );
    }
  }

  return value;
}


double r8mat_norm_fro ( int m, int n, double a[] )












































{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
}


double r8mat_norm_fro_affine ( int m, int n, double a1[], double a2[] )












































{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a1[i+j*m] - a2[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
}


double r8mat_norm_l1 ( int m, int n, double a[] )











































{
  double col_sum;
  int i;
  int j;
  double value;

  value = 0.0;

  for ( j = 0; j < n; j++ )
  {
    col_sum = 0.0;
    for ( i = 0; i < m; i++ )
    {
      col_sum = col_sum + fabs ( a[i+j*m] );
    }
    value = r8_max ( value, col_sum );
  }
  return value;
}


double r8mat_norm_l2 ( int m, int n, double a[] )












































{
  double *at;
  double *b;
  double *diag;
  double value;

  at = r8mat_transpose_new ( m, n, a );



  b = r8mat_mm_new ( m, n, m, a, at );



  r8mat_symm_jacobi ( m, b );



  diag = r8mat_diag_get_vector_new ( m, b );

  value = sqrt ( r8vec_max ( m, diag ) );

  delete [] at;
  delete [] b;
  delete [] diag;

  return value;
}


double r8mat_norm_li ( int m, int n, double a[] )











































{
  int i;
  int j;
  double row_sum;
  double value;

  value = 0.0;

  for ( i = 0; i < m; i++ )
  {
    row_sum = 0.0;
    for ( j = 0; j < n; j++ )
    {
      row_sum = row_sum + fabs ( a[i+j*m] );
    }
    value = r8_max ( value, row_sum );
  }
  return value;
}


double *r8mat_normal_01_new ( int m, int n, int &seed )














































{
  double *r;

  r = r8vec_normal_01_new ( m * n, seed );

  return r;
}


double *r8mat_nullspace ( int m, int n, double a[], int nullspace_size )
























































{
  int *col;
  double det;
  int i;
  int i2;
  int j;
  int j2;
  double *nullspace;
  int *row;
  double *rref;



  rref = r8mat_copy_new ( m, n, a );



  det = r8mat_rref ( m, n, rref );




  row = new int[m];
  for ( i = 0; i < m; i++ )
  {
    row[i] = 0;
  }

  col = new int[n];
  for ( j = 0; j < n; j++ )
  {
    col[j] = - ( j + 1 );
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( rref[i+j*m] == 1.0 )
      {
        row[i] = ( j + 1 );
        col[j] = ( j + 1 );
        break;
      }
    }
  }

  nullspace = r8mat_zeros_new ( n, nullspace_size );

  j2 = 0;




  for ( j = 0; j < n; j++ )
  {
    if ( col[j] < 0 )
    {
      for ( i = 0; i < m; i++ )
      {
        if ( rref[i+j*m] != 0.0 )
        {
          i2 = row[i] - 1;
          nullspace[i2+j2*n] = - rref[i+j*m];
        }
      }
      nullspace[j+j2*n] = 1.0;
      j2 = j2 + 1;
    }
  }
  delete [] col;
  delete [] row;
  delete [] rref;

  return nullspace;
}


int r8mat_nullspace_size ( int m, int n, double a[] )

























































{
  double det;
  int i;
  int j;
  int leading;
  int nullspace_size;
  double *rref;



  rref = r8mat_copy_new ( m, n, a );



  det = r8mat_rref ( m, n, rref );



  leading = 0;
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( rref[i+j*m] == 1.0 )
      {
        leading = leading + 1;
        break;
      }
    }
  }
  nullspace_size = n - leading;

  delete [] rref;

  return nullspace_size;
}


double *r8mat_orth_uniform_new ( int n, int &seed )


















































































{
  double *a_col;
  double *q;
  double *q2;
  int i;
  int j;
  double *v;



  q = r8mat_identity_new ( n );















  a_col = new double[n];

  for ( j = 1; j < n; j++ )
  {



    for ( i = 1; i < j; i++ )
    {
      a_col[i-1] = 0.0;
    }
    for ( i = j; i <= n; i++ )
    {
      a_col[i-1] = r8_normal_01 ( seed );
    }




    v = r8vec_house_column ( n, a_col, j );



    q2 = r8mat_house_axh_new ( n, q, v );

    delete [] v;

    r8mat_copy ( n, n, q2, q );

    delete [] q2;
  }



  delete [] a_col;

  return q;
}


void r8mat_plot ( int m, int n, double a[], string title )


































{
  int i;
  int j;
  int jhi;
  int jlo;

  cout << "\n";
  cout << title << "\n";

  for ( jlo = 1; jlo <= n; jlo = jlo + 70 )
  {
    jhi = i4_min ( jlo + 70-1, n );
    cout << "\n";
    cout << "          ";
    for ( j = jlo; j <= jhi; j++ )
    {
      cout <<  ( j % 10 );
    }
    cout << "\n";
    cout << "\n";

    for ( i = 1; i <= m; i++ )
    {
      cout << setw(6) << i << "    ";
      for ( j = jlo; j <= jhi; j++ )
      {
        cout << r8mat_plot_symbol ( a[i-1+(j-1)*m] );
      }
      cout << "\n";
    }
  }

  return;
}


char r8mat_plot_symbol ( double r )

































{
  char c;

  if ( r < 0.0 )
  {
    c = '-';
  }
  else if ( r == 0.0 )
  {
    c = '0';
  }
  else
  {
    c = '+';
  }

  return c;
}


double *r8mat_poly_char ( int n, double a[] )



































{
  int i;
  int order;
  double *p;
  double trace;
  double *work1;
  double *work2;

  p = new double[n+1];



  work1 = r8mat_identity_new ( n );

  p[n] = 1.0;

  for ( order = n-1; 0 <= order; order-- )
  {



    work2 = r8mat_mm_new ( n, n, n, a, work1 );



    trace = r8mat_trace ( n, work2 );



    p[order] = -trace / ( double ) ( n - order );



    delete [] work1;

    r8mat_copy ( n, n, work2, work1 );

    delete [] work2;

    for ( i = 0; i < n; i++ )
    {
      work1[i+i*n] = work1[i+i*n] + p[order];
    }
  }

  delete [] work1;

  return p;
}


double *r8mat_power ( int n, double a[], int npow )










































{
  double *b;
  double *c;
  int ipow;

  if ( npow < 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_POWER - Fatal error!\n";
    cerr << "  Input value of NPOW < 0.\n";
    cerr << "  NPOW = " << npow << "\n";
    exit ( 1 );
  }

  b = r8mat_identity_new ( n );

  for ( ipow = 1; ipow <= npow; ipow++ )
  {
    c = r8mat_mm_new ( n, n, n, a, b );
    r8mat_copy ( n, n, c, b );
    delete [] c;
  }

  return b;
}


void r8mat_power_method ( int n, double a[], double *r, double v[] )








































{
  double *av;
  double eps;
  int i;
  int it;
  double it_eps = 0.0001;
  int it_max = 100;
  int it_min = 10;
  int j;
  double r2;
  double r_old;

  eps = sqrt ( r8_epsilon ( ) );

  *r = r8vec_norm ( n, v );

  if ( *r == 0.0 )
  {
    for ( i = 0; i < n; i++ )
    {
      v[i] = 1.0;
    }
    *r = sqrt ( ( double ) n );
  }

  for ( i = 0; i < n; i++ )
  {
    v[i] = v[i] / *r;
  }

  for ( it = 1; it <= it_max; it++ )
  {
    av = r8mat_mv_new ( n, n, a, v );

    r_old = *r;
    *r = r8vec_norm ( n, av );

    if ( it_min < it )
    {
      if ( fabs ( *r - r_old ) <= it_eps * ( 1.0 + fabs ( *r ) ) )
      {
        break;
      }
    }

    r8vec_copy ( n, av, v );

    delete [] av;

    if ( *r != 0.0 )
    {
      for ( i = 0; i < n; i++ )
      {
        v[i] = v[i] / *r;
      }
    }




    if ( it < it_max / 2 )
    {
      j = ( ( it - 1 ) % n );
      v[j] = v[j] + eps * ( 1.0 + fabs ( v[j] ) );
      r2 = r8vec_norm ( n, v );
      for ( i = 0; i < n; i++ )
      {
        v[i] = v[i] / r2;
      }
    }
  }
  return;
}


void r8mat_print ( int m, int n, double a[], string title )




































{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}


void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title )







































{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }



  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }
    cout << "\n";





    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j - 1 << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";



    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( ihi < m )
    {
      i2hi = ihi;
    }
    else
    {
      i2hi = m;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {



      cout << setw(5) << i - 1 << ": ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}


double r8mat_product_elementwise ( int m, int n, double a[], double b[] )



































{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + a[i+j*m] * b[i+j*m];
    }
  }
  
  return value;
}


double r8mat_ref ( int m, int n, double a[] )



























































{
  double asum;
  double det;
  int i;
  int j;
  int lead;
  int r;
  double temp;
  double tol;

  det = 1.0;
  asum = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      asum = asum + fabs ( a[i+j*m] );
    }
  }
  tol = r8_epsilon ( ) * asum;
  lead = 0;

  for ( r = 0; r < m; r++ )
  {
    if ( n - 1 < lead )
    {
      break;
    }

    i = r;

    while ( fabs ( a[i+lead*m] ) <= tol )
    {
      i = i + 1;

      if ( m - 1 < i )
      {
        i = r;
        lead = lead + 1;
        if ( n - 1 < lead )
        {
          lead = -1;
          break;
         }
      }
    }

    if ( lead < 0 )
    {
      break;
    }

    for ( j = 0; j < n; j++ )
    {
      temp     = a[i+j*m];
      a[i+j*m] = a[r+j*m];
      a[r+j*m] = temp;
    }

    det = det * a[r+lead*m];
    temp = a[r+lead*m];

    for ( j = 0; j < n; j++ )
    {
      a[r+j*m] = a[r+j*m] / temp;
    }

    for ( i = r + 1; i < m; i++ )
    {
      temp = a[i+lead*m];
      for ( j = 0; j < n; j++ )
      {
        a[i+j*m] = a[i+j*m] - temp * a[r+j*m];
      }
    }
    lead = lead + 1;
  }
  return det;
}


double r8mat_rms ( int m, int n, double a[] )




































{
  int i;
  int j;
  double value;

  value = 0.0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + a[i+j*m] * a[i+j*m];
    }
    value = sqrt ( value / ( double ) ( m ) / ( double ) ( n ) );
  }
  return value;
}


void r8mat_row_copy ( int m, int n, int i, double v[], double a[] )



































{
  int j;

  for ( j = 0; j < n; j++ )
  {
    a[i+j*m] = v[j];
  }
  return;
}


double r8mat_rref ( int m, int n, double a[] )
































































{
  double asum;
  double det;
  int i;
  int j;
  int lead;
  int r;
  double temp;
  double tol;

  det = 1.0;
  asum = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      asum = asum + fabs ( a[i+j*m] );
    }
  }
  tol = r8_epsilon ( ) * asum;
  lead = 0;

  for ( r = 0; r < m; r++ )
  {
    if ( n - 1 < lead )
    {
      break;
    }

    i = r;

    while ( fabs ( a[i+lead*m] ) <= tol )
    {
      i = i + 1;

      if ( m - 1 < i )
      {
        i = r;
        lead = lead + 1;
        if ( n - 1 < lead )
        {
          lead = -1;
          break;
         }
      }
    }

    if ( lead < 0 )
    {
      break;
    }

    for ( j = 0; j < n; j++ )
    {
      temp     = a[i+j*m];
      a[i+j*m] = a[r+j*m];
      a[r+j*m] = temp;
    }

    det = det * a[r+lead*m];
    temp = a[r+lead*m];

    for ( j = 0; j < n; j++ )
    {
      a[r+j*m] = a[r+j*m] / temp;
    }

    for ( i = 0; i < m; i++ )
    {
      if ( i != r )
      {
        temp = a[i+lead*m];
        for ( j = 0; j < n; j++ )
        {
          a[i+j*m] = a[i+j*m] - temp * a[r+j*m];
        }
      }
    }
    lead = lead + 1;

  }
  return det;
}


void r8mat_scale ( int m, int n, double s, double a[] )































{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] * s;
    }
  }
  return;
}


int r8mat_solve ( int n, int rhs_num, double a[] )












































{
  double apivot;
  double factor;
  int i;
  int ipivot;
  int j;
  int k;
  double temp;

  for ( j = 0; j < n; j++ )
  {



    ipivot = j;
    apivot = a[j+j*n];

    for ( i = j; i < n; i++ )
    {
      if ( fabs ( apivot ) < fabs ( a[i+j*n] ) )
      {
        apivot = a[i+j*n];
        ipivot = i;
      }
    }

    if ( apivot == 0.0 )
    {
      return j;
    }



    for ( i = 0; i < n + rhs_num; i++ )
    {
      temp          = a[ipivot+i*n];
      a[ipivot+i*n] = a[j+i*n];
      a[j+i*n]      = temp;
    }



    a[j+j*n] = 1.0;
    for ( k = j; k < n + rhs_num; k++ )
    {
      a[j+k*n] = a[j+k*n] / apivot;
    }



    for ( i = 0; i < n; i++ )
    {
      if ( i != j )
      {
        factor = a[i+j*n];
        a[i+j*n] = 0.0;
        for ( k = j; k < n + rhs_num; k++ )
        {
          a[i+k*n] = a[i+k*n] - factor * a[j+k*n];
        }
      }
    }
  }

  return 0;
}


double *r8mat_solve_2d ( double a[], double b[], double *det )










































{
  double *x;



  *det = a[0+0*2] * a[1+1*2] - a[0+1*2] * a[1+0*2];



  if ( *det == 0.0 )
  {
    return NULL;
  }



  x = new double[2];

  x[0] = (  a[1+1*2] * b[0] - a[0+1*2] * b[1] ) / ( *det );
  x[1] = ( -a[1+0*2] * b[0] + a[0+0*2] * b[1] ) / ( *det );

  return x;
}


double *r8mat_solve_3d ( double a[], double b[], double *det )










































{
  double *x;



  *det =  a[0+0*3] * ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] )
        + a[0+1*3] * ( a[1+2*3] * a[2+0*3] - a[1+0*3] * a[2+2*3] )
        + a[0+2*3] * ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] );



  if ( *det == 0.0 )
  {
    return NULL;
  }



  x = new double[3];

  x[0] = (   ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] ) * b[0]
           - ( a[0+1*3] * a[2+2*3] - a[0+2*3] * a[2+1*3] ) * b[1]
           + ( a[0+1*3] * a[1+2*3] - a[0+2*3] * a[1+1*3] ) * b[2] ) / ( *det );

  x[1] = ( - ( a[1+0*3] * a[2+2*3] - a[1+2*3] * a[2+0*3] ) * b[0]
           + ( a[0+0*3] * a[2+2*3] - a[0+2*3] * a[2+0*3] ) * b[1]
           - ( a[0+0*3] * a[1+2*3] - a[0+2*3] * a[1+0*3] ) * b[2] ) / ( *det );

  x[2] = (   ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] ) * b[0]
           - ( a[0+0*3] * a[2+1*3] - a[0+1*3] * a[2+0*3] ) * b[1]
           + ( a[0+0*3] * a[1+1*3] - a[0+1*3] * a[1+0*3] ) * b[2] ) / ( *det );

  return x;
}


double *r8mat_solve2 ( int n, double a[], double b[], int &ierror )


















































{
  double amax;
  int i;
  int imax;
  int j;
  int k;
  int *piv;
  double *x;

  ierror = 0;

  piv = i4vec_zeros_new ( n );
  x = r8vec_zeros_new ( n );



  for ( k = 1; k <= n; k++ )
  {






    amax = 0.0;
    imax = 0;
    for ( i = 1; i <= n; i++ )
    {
      if ( piv[i-1] == 0 )
      {
        if ( amax < fabs ( a[i-1+(k-1)*n] ) )
        {
          imax = i;
          amax = fabs ( a[i-1+(k-1)*n] );
        }
      }
    }




    if ( imax != 0 )
    {
      piv[imax-1] = k;
      for ( j = k+1; j <= n; j++ )
      {
        a[imax-1+(j-1)*n] = a[imax-1+(j-1)*n] / a[imax-1+(k-1)*n];
      }
      b[imax-1] = b[imax-1] / a[imax-1+(k-1)*n];
      a[imax-1+(k-1)*n] = 1.0;

      for ( i = 1; i <= n; i++ )
      {
        if ( piv[i-1] == 0 )
        {
          for ( j = k+1; j <= n; j++ )
          {
            a[i-1+(j-1)*n] = a[i-1+(j-1)*n] - a[i-1+(k-1)*n] * a[imax-1+(j-1)*n];
          }
          b[i-1] = b[i-1] - a[i-1+(k-1)*n] * b[imax-1];
          a[i-1+(k-1)*n] = 0.0;
        }
      }
    }
  }




  for ( j = n; 1 <= j; j-- )
  {
    imax = 0;
    for ( k = 1; k <= n; k++ )
    {
      if ( piv[k-1] == j )
      {
        imax = k;
      }
    }

    if ( imax == 0 )
    {
      x[j-1] = 0.0;

      if ( b[j-1] == 0.0 )
      {
        ierror = 1;
        cout << "\n";
        cout << "R8MAT_SOLVE2 - Warning:\n";
        cout << "  Consistent singularity, equation = " << j << "\n";
      }
      else
      {
        ierror = 2;
        cout << "\n";
        cout << "R8MAT_SOLVE2 - Warning:\n";
        cout << "  Inconsistent singularity, equation = " << j << "\n";
      }
    }
    else
    {
      x[j-1] = b[imax-1];

      for ( i = 1; i <= n; i++ )
      {
        if ( i != imax )
        {
          b[i-1] = b[i-1] - a[i-1+(j-1)*n] * x[j-1];
        }
      }
    }
  }

  delete [] piv;

  return x;
}


double r8mat_sum ( int m, int n, double a[] )
































{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + a[i+j*m];
    }
  }
  return value;
}


double *r8mat_symm_eigen ( int n, double x[], double q[] )








































{
  double *a;
  int i;
  int j;
  int k;



  a = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i+j*n] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        a[i+j*n] = a[i+j*n] + q[i+k*n] * x[k] * q[j+k*n];
      }
    }
  }

  return a;
}


void r8mat_symm_jacobi ( int n, double a[] )



































{
  double c;
  double eps = 0.00001;
  int i;
  int it;
  int it_max = 100;
  int j;
  int k;
  double norm_fro;
  double s;
  double sum2;
  double t;
  double t1;
  double t2;
  double u;

  norm_fro = r8mat_norm_fro ( n, n, a );

  it = 0;

  for ( ; ; )
  {
    it = it + 1;

    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < i; j++ )
      {
        if ( eps * norm_fro < fabs ( a[i+j*n] ) + fabs ( a[j+i*n] ) )
        {
          u = ( a[j+j*n] - a[i+i*n] ) / ( a[i+j*n] + a[j+i*n] );

          t = r8_sign ( u ) / ( fabs ( u ) + sqrt ( u * u + 1.0 ) );
          c = 1.0 / sqrt ( t * t + 1.0 );
          s = t * c;



          for ( k = 0; k < n; k++ )
          {
            t1 = a[i+k*n];
            t2 = a[j+k*n];
            a[i+k*n] = t1 * c - t2 * s;
            a[j+k*n] = t1 * s + t2 * c;
          }



          for ( k = 0; k < n; k++ )
          {
            t1 = a[k+i*n];
            t2 = a[k+j*n];
            a[k+i*n] = c * t1 - s * t2;
            a[k+j*n] = s * t1 + c * t2;
          }
        }
      }
    }



    sum2 = 0.0;
    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < i; j++ )
      {
        sum2 = sum2 + fabs ( a[i+j*n] );
      }
    }

    if ( sum2 <= eps * ( norm_fro + 1.0 ) )
    {
      break;
    }

    if ( it_max <= it )
    {
      break;
    }

  }

  return;
}


double **r8mat_to_r8cmat_new (  int m, int n, double a[] )






































{
  double **b;
  int i;
  int j;

  b = r8cmat_new ( m, n );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      b[j][i] = a[i+j*m];
    }
  }

  return b;
}


int r8mat_to_r8plu ( int n, double a[], int pivot[], double lu[] )



















































{
  int i;
  int info;
  int j;
  int k;
  int l;
  double temp;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      lu[i+j*n] = a[i+j*n];
    }
  }
  info = 0;

  for ( k = 1; k <= n-1; k++ )
  {



    l = k;
    for ( i = k+1; i <= n; i++ )
    {
      if ( fabs ( lu[l-1+(k-1)*n] ) < fabs ( lu[i-1+(k-1)*n] ) )
      {
        l = i;
      }
    }

    pivot[k-1] = l;



    if ( lu[l-1+(k-1)*n] == 0.0 )
    {
      info = k;
      return info;
    }



    if ( l != k )
    {
      temp            = lu[l-1+(k-1)*n];
      lu[l-1+(k-1)*n] = lu[k-1+(k-1)*n];
      lu[k-1+(k-1)*n] = temp;
    }



    for ( i = k+1; i <= n; i++ )
    {
      lu[i-1+(k-1)*n] = -lu[i-1+(k-1)*n] / lu[k-1+(k-1)*n];
    }



    for ( j = k+1; j <= n; j++ )
    {
      if ( l != k )
      {
        temp            = lu[l-1+(j-1)*n];
        lu[l-1+(j-1)*n] = lu[k-1+(j-1)*n];
        lu[k-1+(j-1)*n] = temp;
      }

      for ( i = k+1; i <= n; i++ )
      {
        lu[i-1+(j-1)*n] = lu[i-1+(j-1)*n] + lu[i-1+(k-1)*n] * lu[k-1+(j-1)*n];
      }
    }
  }

  pivot[n-1] = n;

  if ( lu[n-1+(n-1)*n] == 0.0 )
  {
    info = n;
  }

  return info;
}


double **r8mat_to_r8rmat ( int m, int n, double a[] )







































{
  double **b;
  int i;
  int j;

  b = r8rmat_new ( m, n );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      b[i][j] = a[i+j*m];
    }
  }

  return b;
}


double r8mat_trace ( int n, double a[] )


































{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i+i*n];
  }

  return value;
}


void r8mat_transpose_in_place ( int n, double a[] )






























{
  int i;
  int j;
  double t;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      t        = a[i+j*n];
      a[i+j*n] = a[j+i*n];
      a[j+i*n] = t;
    }
  }
  return;
}


double *r8mat_transpose_new ( int m, int n, double a[] )
































{
  double *b;
  int i;
  int j;

  b = new double[n*m];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      b[j+i*n] = a[i+j*m];
    }
  }
  return b;
}


void r8mat_transpose_print ( int m, int n, double a[], string title )
































{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}


void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title )




































{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int i2lo_hi;
  int i2lo_lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }

  if ( ilo < 1 )
  {
    i2lo_lo = 1;
  }
  else
  {
    i2lo_lo = ilo;
  }

  if ( ihi < m )
  {
    i2lo_hi = m;
  }
  else
  {
    i2lo_hi = ihi;
  }

  for ( i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;

    if ( m < i2hi )
    {
      i2hi = m;
    }
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i - 1 << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    if ( jlo < 1 )
    {
      j2lo = 1;
    }
    else
    {
      j2lo = jlo;
    }
    if ( n < jhi )
    {
      j2hi = n;
    }
    else
    {
      j2hi = jhi;
    }

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j - 1 << ":";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}


double *r8mat_u_inverse ( int n, double a[] )













































{
  double *b;
  int i;
  int j;
  int k;

  b = new double[n*n];

  for ( j = n-1; 0 <= j; j-- )
  {
    for ( i = n-1; 0 <= i; i-- )
    {
      if ( j < i )
      {
        b[i+j*n] = 0.0;
      }
      else if ( i == j )
      {
        b[i+j*n] = 1.0 / a[i+j*n];
      }
      else
      {
        b[i+j*n] = 0.0;
        for ( k = i+1; k <= j; k++ )
        {
          b[i+j*n] = b[i+j*n] - a[i+k*n] * b[k+j*n];
        }
       b[i+j*n] = b[i+j*n] / a[i+i*n];
      }
    }
  }

  return b;
}


double *r8mat_u_solve ( int n, double a[], double b[] )


































{
  int i;
  int j;
  double *x;



  x = new double[n];

  for ( i = n - 1; 0 <= i; i-- )
  {
    x[i] = b[i];
    for ( j = i + 1; j < n; j++ )
    {
      x[i] = x[i] - a[i+j*n] * x[j];
    }
    x[i] = x[i] / a[i+i*n];
  }

  return x;
}


double *r8mat_u1_inverse ( int n, double a[] )













































{
  double *b;
  int i;
  int j;
  int k;

  b = new double[n*n];

  for ( j = n-1; 0 <= j; j-- )
  {
    for ( i = n-1; 0 <= i; i-- )
    {
      if ( j < i )
      {
        b[i+j*n] = 0.0;
      }
      else if ( i == j )
      {
        b[i+j*n] = 1.0;
      }
      else
      {
        b[i+j*n] = 0.0;
        for ( k = i+1; k <= j; k++ )
        {
          b[i+j*n] = b[i+j*n] - a[i+k*n] * b[k+j*n];
        }
       b[i+j*n] = b[i+j*n] / a[i+i*n];
      }
    }
  }

  return b;
}


void r8mat_uniform_01 ( int m, int n, int &seed, double r[] )






































































{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }

      r[i+j*m] = ( double ) ( seed ) * 4.656612875E-10;
    }
  }
  return;
}


double *r8mat_uniform_01_new ( int m, int n, int &seed )





























































{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;
  double *r;

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }
      r[i+j*m] = ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return r;
}


void r8mat_uniform_ab ( int m, int n, double a, double b, int &seed, double r[] )








































































{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }

      r[i+j*m] = a + ( b - a ) * ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return;
}


double *r8mat_uniform_ab_new ( int m, int n, double a, double b, int &seed )








































































{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_AB_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }

      r[i+j*m] = a + ( b - a ) * ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return r;
}


void r8mat_uniform_abvec ( int m, int n, double a[], double b[], int &seed, 
  double r[] )








































































{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_ABVEC - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }

      r[i+j*m] = a[i] + ( b[i] - a[i] ) * ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return;
}


double *r8mat_uniform_abvec_new ( int m, int n, double a[], double b[], 
  int &seed )









































































{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_ABVEC_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }

      r[i+j*m] = a[i] + ( b[i] - a[i] ) * ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return r;
}


double *r8mat_ut_solve ( int n, double a[], double b[] )






































{
  int i;
  int j;
  double *x;



  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
    for ( j = 0; j < i; j++ )
    {
      x[i] = x[i] - a[j+i*n] * x[j];
    }
    x[i] = x[i] / a[i+i*n];
  }

  return x;
}


double *r8mat_vand2 ( int n, double x[] )





























































{
  double *a;
  int i;
  int j;

  a = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( j == 0 && x[i] == 0.0 )
      {
        a[i+j*n] = 1.0;
      }
      else
      {
        a[i+j*n] = pow ( x[i], j );
      }
    }
  }

  return a;
}


double r8mat_vtmv ( int m, int n, double x[], double a[], double y[] )




































{
  int i;
  int j;
  double vtmv;

  vtmv = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      vtmv = vtmv + x[i] * a[i+j*m] * y[j];
    }
  }
  return vtmv;
}


void r8mat_zeros ( int m, int n, double a[] )






























{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  return;
}


double *r8mat_zeros_new ( int m, int n )






























{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  return a;
}


double r8plu_det ( int n, int pivot[], double lu[] )








































{
  double det;
  int i;

  det = 1.0;

  for ( i = 0; i < n; i++ )
  {
    det = det * lu[i+i*n];
    if ( pivot[i] != i+1 )
    {
      det = -det;
    }
  }

  return det;
}


void r8plu_inverse ( int n, int pivot[], double lu[], double a_inverse[] )


































{
  int i;
  int j;
  int k;
  double temp;
  double *work;

  work = new double[n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a_inverse[i+j*n] = lu[i+j*n];
    }
  }



  for ( k = 1; k <= n; k++ )
  {
    a_inverse[k-1+(k-1)*n]     = 1.0 / a_inverse[k-1+(k-1)*n];
    for ( i = 1; i <= k-1; i++ )
    {
      a_inverse[i-1+(k-1)*n] = -a_inverse[i-1+(k-1)*n] * a_inverse[k-1+(k-1)*n];
    }

    for ( j = k+1; j <= n; j++ )
    {
      temp                     = a_inverse[k-1+(j-1)*n];
      a_inverse[k-1+(j-1)*n]   = 0.0;
      for ( i = 1; i <= k; i++ )
      {
        a_inverse[i-1+(j-1)*n] = a_inverse[i-1+(j-1)*n]
          + temp * a_inverse[i-1+(k-1)*n];
      }
    }
  }



  for ( k = n-1; 1 <= k; k-- )
  {
    for ( i = k+1; i <= n; i++ )
    {
      work[i-1] = a_inverse[i-1+(k-1)*n];
      a_inverse[i-1+(k-1)*n] = 0.0;
    }

    for ( j = k+1; j <= n; j++ )
    {
      for ( i = 1; i <= n; i++ )
      {
        a_inverse[i-1+(k-1)*n] = a_inverse[i-1+(k-1)*n]
          + a_inverse[i-1+(j-1)*n] * work[j-1];
      }
    }

    if ( pivot[k-1] != k )
    {
      for ( i = 1; i <= n; i++ )
      {
        temp                            = a_inverse[i-1+(k-1)*n];
        a_inverse[i-1+(k-1)*n]          = a_inverse[i-1+(pivot[k-1]-1)*n];
        a_inverse[i-1+(pivot[k-1]-1)*n] = temp;
      }
    }
  }

  delete [] work;

  return;
}


void r8plu_mul ( int n, int pivot[], double lu[], double x[], double b[] )





































{
  int i;
  int j;
  int k;
  double temp;

  for ( i = 0; i < n; i++ )
  {
    b[i] = x[i];
  }



  for ( j = 1; j <= n; j++ )
  {
    for ( i = 0; i < j-1; i++ )
    {
      b[i] = b[i] + lu[i+(j-1)*n] * b[j-1];
    }
    b[j-1] = lu[j-1+(j-1)*n] * b[j-1];
  }



  for ( j = n-1; 1 <= j; j-- )
  {
    for ( i = j; i < n; i++ )
    {
      b[i] = b[i] - lu[i+(j-1)*n] * b[j-1];
    }

    k = pivot[j-1];

    if ( k != j )
    {
      temp = b[k-1];
      b[k-1] = b[j-1];
      b[j-1] = temp;
    }
  }

  return;
}


void r8plu_sol ( int n, int pivot[], double lu[], double b[], double x[] )



































{
  int i;
  int j;
  int k;
  double temp;



  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  for ( k = 1; k <= n-1; k++ )
  {
    j = pivot[k-1];

    if ( j != k )
    {
      temp   = x[j-1];
      x[j-1] = x[k-1];
      x[k-1] = temp;
    }

    for ( i = k+1; i <= n; i++ )
    {
      x[i-1] = x[i-1] + lu[i-1+(k-1)*n] * x[k-1];
    }
  }



  for ( k = n; 1 <= k; k-- )
  {
    x[k-1] = x[k-1] / lu[k-1+(k-1)*n];
    for ( i = 1; i <= k-1; i++ )
    {
      x[i-1] = x[i-1] - lu[i-1+(k-1)*n] * x[k-1];
    }
  }

  return;
}


void r8plu_to_r8mat ( int n, int pivot[], double lu[], double a[] )































{
  int i;
  int j;
  int k;
  double temp;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = 1.0;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }

  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= n; i++ )
    {
      for ( k = 1; k <= i-1; k++ )
      {
        a[k-1+(j-1)*n] = a[k-1+(j-1)*n] + lu[k-1+(i-1)*n] * a[i-1+(j-1)*n];
      }
      a[i-1+(j-1)*n] = lu[i-1+(i-1)*n] * a[i-1+(j-1)*n];
    }



    for ( i = n-1; 1 <= i; i-- )
    {
      for ( k = i+1; k <= n; k++ )
      {
        a[k-1+(j-1)*n] = a[k-1+(j-1)*n] - lu[k-1+(i-1)*n] * a[i-1+(j-1)*n];
      }

      k = pivot[i-1];

      if ( k != i )
      {
        temp           = a[k-1+(j-1)*n];
        a[k-1+(j-1)*n] = a[i-1+(j-1)*n];
        a[i-1+(j-1)*n] = temp;
      }
    }
  }

  return;
}


int r8poly_degree ( int na, double a[] )




































{
  int degree;

  degree = na;

  while ( 0 < degree )
  {
    if ( a[degree] != 0.0 )
    {
      return degree;
    }
    degree = degree - 1;
  }

  return degree;
}


double *r8poly_deriv ( int n, double c[], int p )





































{
  double *cp;
  double *cp_temp;
  int d;
  int i;

  if ( n < p )
  {
    return NULL;
  }
  cp_temp = r8vec_copy_new ( n+1, c );

  for ( d = 1; d <= p; d++ )
  {
    for ( i = 0; i <= n-d; i++ )
    {
      cp_temp[i] = ( double ) ( i + 1 ) * cp_temp[i+1];
    }
    cp_temp[n-d+1] = 0.0;
  }

  cp = r8vec_copy_new ( n - p + 1, cp_temp );

  delete [] cp_temp;

  return cp;
}


double r8poly_lagrange_0 ( int npol, double xpol[], double xval )














































{
  int i;
  double wval;

  wval = 1.0;
  for ( i = 0; i < npol; i++ )
  {
    wval = wval * ( xval - xpol[i] );
  }

  return wval;
}


double r8poly_lagrange_1 ( int npol, double xpol[], double xval )











































{
  double dwdx;
  int i;
  double w;

  dwdx = 0.0;
  w = 1.0;

  for ( i = 0; i < npol; i++ )
  {
    dwdx = w + ( xval - xpol[i] ) * dwdx;
    w = w * ( xval - xpol[i] );
  }

  return dwdx;
}


double r8poly_lagrange_2 ( int npol, double xpol[], double xval )



















































{
  double dw2dx2;
  int i;
  int j;
  int k;
  double term;

  dw2dx2 = 0.0;

  for ( k = 0; k < npol; k++ )
  {
    for ( j = 0; j < npol; j++ )
    {
      if ( j != k )
      {
        term = 1.0;
        for ( i = 0; i < npol; i++ )
        {
          if ( i != j && i != k )
          {
            term = term * ( xval - xpol[i] );
          }
        }
        dw2dx2 = dw2dx2 + term;
      }
    }
  }

  return dw2dx2;
}


double *r8poly_lagrange_coef ( int npol, int ipol, double xpol[] )








































{
  int i;
  int index;
  int j;
  double *pcof;



  if ( ipol < 1 || npol < ipol )
  {
    cerr << "\n";
    cerr << "R8POLY_LAGRANGE_COEF - Fatal error!\n";
    cerr << "  1 <= IPOL <= NPOL is required.\n";
    cerr << "  but IPOL = " << ipol << "\n";
    cerr << "  and NPOL = " << npol << "\n";
    exit ( 1 );
  }



  if ( ! r8vec_is_distinct ( npol, xpol ) )
  {
    cerr << "\n";
    cerr << "R8POLY_LAGRANGE_COEF - Fatal error!\n";
    cerr << "  Two entries of XPOL are equal:\n";
    exit ( 1 );
  }

  pcof = new double[npol];

  pcof[0] = 1.0;
  for ( i = 1; i < npol; i++ )
  {
    pcof[i] = 0.0;
  }

  index = 0;

  for ( i = 1; i <= npol; i++ )
  {
    if ( i != ipol )
    {
      index = index + 1;

      for ( j = index; 0 <= j; j-- )
      {
        pcof[j] = - xpol[i-1] * pcof[j] / ( xpol[ipol-1] - xpol[i-1] );

        if ( 0 < j )
        {
          pcof[j] = pcof[j] + pcof[j-1] / ( xpol[ipol-1] - xpol[i-1] );
        }
      }
    }
  }

  return pcof;
}


void r8poly_lagrange_factor ( int npol, double xpol[], double xval,
  double *wval, double *dwdx )









































































{
  int i;
  int j;
  double term;

  *wval = 1.0;
  for ( i = 0; i < npol; i++ )
  {
    *wval = *wval * ( xval - xpol[i] );
  }

  *dwdx = 0.0;

  for ( i = 0; i < npol; i++ )
  {
    term = 1.0;

    for ( j = 0; j < npol; j++ )
    {
      if ( i != j )
      {
        term = term * ( xval - xpol[j] );
      }
    }
    *dwdx = *dwdx + term;
  }

  return;
}


int r8poly_lagrange_val ( int npol, int ipol, double xpol[], double xval,
  double *pval, double *dpdx )














































{
  int i;
  int j;
  double p2;



  if ( ipol < 0 || npol-1 < ipol )
  {
    cerr << "\n";
    cerr << "R8POLY_LAGRANGE_VAL - Fatal error!\n";
    cerr << "  0 <= IPOL <= NPOL-1 is required.\n";
    exit ( 1 );
  }



  for ( i = 1; i < npol; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      if ( xpol[i] == xpol[j] )
      {
        cerr << "\n";
        cerr << "R8POLY_LAGRANGE_VAL - Fatal error!\n";
        cerr << "  Two entries of XPOL are equal:\n";
        cerr << "  XPOL(" << i << ") = " << xpol[i] << ".\n";
        cerr << "  XPOL(" << j << ") = " << xpol[j] << ".\n";
        exit ( 1 );
      }
    }
  }



  *pval = 1.0;

  for ( i = 0; i < npol; i++ )
  {
    if ( i != ipol )
    {
      *pval = *pval * ( xval - xpol[i] ) / ( xpol[ipol] - xpol[i] );
    }
  }




  *dpdx = 0.0;

  for ( i = 0; i < npol; i++ )
  {
    if ( i != ipol )
    {
      p2 = 1.0;

      for ( j = 0; j < npol; j++ )
      {
        if ( j == i )
        {
          p2 = p2                      / ( xpol[ipol] - xpol[j] );
        }
        else if ( j != ipol )
        {
          p2 = p2 * ( xval - xpol[j] ) / ( xpol[ipol] - xpol[j] );
        }
      }
      *dpdx = *dpdx + p2;
    }
  }

  return 0;
}


int r8poly_order ( int na, double a[] )



































{
  int order;

  order = na + 1;

  while ( 1 < order )
  {
    if ( a[order-1] != 0.0 )
    {
      return order;
    }
    order = order - 1;
  }

  return order;
}


void r8poly_print ( int n, double a[], string title )





























{
  int i;
  double mag;
  char plus_minus;

  if ( 0 < title.length ( ) )
  {
    cout << "\n";
    cout << title << "\n";
  }
  cout << "\n";

  if ( n < 0 )
  {
    cout << "  p(x) = 0\n";
    return;
  }

  if ( a[n] < 0.0 )
  {
    plus_minus = '-';
  }
  else
  {
    plus_minus = ' ';
  }

  mag = fabs ( a[n] );

  if ( 2 <= n )
  {
    cout << "  p(x) = " << plus_minus
         << setw(14) << mag << " * x ^ " << n << "\n";
  }
  else if ( n == 1 )
  {
    cout << "  p(x) = " << plus_minus
         << setw(14) << mag << " * x\n";
  }
  else if ( n == 0 )
  {
    cout << "  p(x) = " << plus_minus
         << setw(14) << mag << "\n";
  }

  for ( i = n - 1; 0 <= i; i-- )
  {
    if ( a[i] < 0.0 )
    {
      plus_minus = '-';
    }
    else
    {
      plus_minus = '+';
    }

    mag = fabs ( a[i] );

    if ( mag != 0.0 )
    {
      if ( 2 <= i )
      {
        cout << "         " << plus_minus
             << setw(14) << mag << " * x ^ " << i << "\n";
      }
      else if ( i == 1 )
      {
        cout << "         " << plus_minus
             << setw(14) << mag << " * x\n";
      }
      else if ( i == 0 )
      {
        cout << "         " << plus_minus
             << setw(14) << mag << "\n";
      }
    }
  }

  return;
}


void r8poly_shift ( double scale, double shift, int n, double poly_cof[] )










































































{
  int i;
  int j;

  for ( i = 1; i <= n; i++ )
  {
    for ( j = i; j <= n; j++ )
    {
      poly_cof[j] = poly_cof[j] / scale;
    }
  }

  for ( i = 0; i <= n - 1; i++ )
  {
    for ( j = n - 1; i <= j; j-- )
    {
      poly_cof[j] = poly_cof[j] - shift * poly_cof[j+1];
    }
  }

  return;
}


double r8poly_value ( int m, double c[], double x )






































{
  int i;
  double value;
  double xi;

  value = c[0];
  xi = 1.0;

  for ( i = 1; i <= m; i++ )
  {
    xi = xi * x;
    value = value + c[i] * xi;
  }

  return value;
}


double r8poly_value_horner ( int m, double c[], double x )






































{
  int i;
  double value;

  value = c[m];

  for ( i = m - 1; 0 <= i; i-- )
  {
    value = value * x + c[i];
  }

  return value;
}


double *r8poly_values_horner ( int m, double c[], int n, double x[] )








































{
  int i;
  int j;
  double *p;

  p = new double[n];

  for ( j = 0; j < n; j++ )
  {
    p[j] = c[m];
  }

  for ( i = m - 1; 0 <= i; i-- )
  {
    for ( j = 0; j < n; j++ )
    {
      p[j] = p[j] * x[j] + c[i];
    }
  }
  return p;
}


double *r8poly_value_2d ( int m, double c[], int n, double x[], double y[] )













































{
  int ex;
  int ey;
  int i;
  int j;
  double *p;
  int s;

  p = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    p[i] = 0.0;
  }

  j = 0;
  for ( s = 0; s <= m; s++ )
  {
    for ( ex = s; 0 <= ex; ex-- )
    {
      ey = s - ex;
      for ( i = 0; i < n; i++ )
      {
        p[i] = p[i] + c[j] * pow ( x[i], ex ) * pow ( y[i], ey );
      }
      j = j + 1;
    }
  }
  return p;
}


int r8poly2_ex ( double x1, double y1, double x2, double y2, double x3,
  double y3, double *x, double *y )


































{
  double bot;

  *x = 0.0;
  *y = 0.0;

  if ( x1 == x2 || x2 == x3 || x3 == x1 )
  {
    return 1;
  }

  if ( y1 == y2 && y2 == y3 && y3 == y1 )
  {
    *x = x1;
    *y = y1;
    return 3;
  }

  bot = ( x2 - x3 ) * y1 + ( x3 - x1 ) * y2 + ( x1 - x2 ) * y3;

  if ( bot == 0.0 )
  {
    return 2;
  }

  *x = 0.5 * (
      x1 * x1 * ( y3 - y2 )
    + x2 * x2 * ( y1 - y3 )
    + x3 * x3 * ( y2 - y1 ) ) /
    ( ( x2 - x3 ) * y1 + ( x3 - x1 ) * y2 + ( x1 - x2 ) * y3 );

  *y = - (
      ( *x - x2 ) * ( *x - x3 ) * ( x2 - x3 ) * y1
    + ( *x - x1 ) * ( *x - x3 ) * ( x3 - x1 ) * y2
    + ( *x - x1 ) * ( *x - x2 ) * ( x1 - x2 ) * y3 ) /
    ( ( x1 - x2 ) * ( x2 - x3 ) * ( x3 - x1 ) );

  return 0;
}


int r8poly2_ex2 ( double x1, double y1, double x2, double y2, double x3,
  double y3, double *x, double *y, double *a, double *b, double *c )





































{
  double v[3*3];
  double *w;

  *a = 0.0;
  *b = 0.0;
  *c = 0.0;
  *x = 0.0;
  *y = 0.0;

  if ( x1 == x2 || x2 == x3 || x3 == x1 )
  {
    return 1;
  }

  if ( y1 == y2 && y2 == y3 && y3 == y1 )
  {
    *x = x1;
    *y = y1;
    return 3;
  }



  v[0+0*3] = 1.0;
  v[0+1*3] = x1;
  v[0+2*3] = x1 * x1;

  v[1+0*3] = 1.0;
  v[1+1*3] = x2;
  v[1+2*3] = x2 * x2;

  v[2+0*3] = 1.0;
  v[2+1*3] = x3;
  v[2+2*3] = x3 * x3;



  w = r8mat_inverse_3d ( v );



  *c = w[0+0*3] * y1 + w[0+1*3] * y2 + w[0+2*3] * y3;
  *b = w[1+0*3] * y1 + w[1+1*3] * y2 + w[1+2*3] * y3;
  *a = w[2+0*3] * y1 + w[2+1*3] * y2 + w[2+2*3] * y3;



  if ( *a == 0.0 )
  {
    return 2;
  }

  *x = - *b / ( 2.0 * *a );
  *y = *a * *x * *x + *b * *x + *c;

  return 0;
}


void r8poly2_rroot ( double a, double b, double c, double *r1, double *r2 )




































{
  double disc;
  double q;

  if ( a == 0.0 )
  {
    cerr << "\n";
    cerr << "R8POLY2_RROOT - Fatal error!\n";
    cerr << "  The coefficient A is zero.\n";
    exit ( 1 );
  }

  disc = b * b - 4.0 * a * c;
  if ( 0.0 <= disc )
  {
    q = ( b + r8_sign ( b ) * sqrt ( disc ) );
    *r1 = -0.5 * q / a;
    *r2 = -2.0 * c / q;
  }
  else
  {
    *r1 = b / 2.0 / a;
    *r2 = b / 2.0 / a;
  }

  return;
}


void r8poly2_val ( double x1, double y1, double x2, double y2,
  double x3, double y3, double x, double *y, double *yp, double *ypp )






































{
  int distinct;
  double dif1;
  double dif2 = 0.0;
  double temp;



  if ( x1 == x2 && x2 == x3 )
  {
    distinct = 1;
  }
  else if ( x1 == x2 )
  {
    distinct = 2;
  }
  else if ( x1 == x3 )
  {
    cerr << "\n";
    cerr << "R8POLY2_VAL - Fatal error!\n";
    cerr << "  X1 = X3 =/= X2.\n";
    return;
  }
  else if ( x2 == x3 )
  {
    distinct = 2;
    temp = x1;
    x1 = x3;
    x3 = temp;
    temp = y1;
    y1 = y2;
    y2 = y3;
    y3 = y1;
  }
  else
  {
    distinct = 3;
  }



  if ( distinct == 1 )
  {
    dif1 = y2;
    dif2 = 0.5 * y3;
  }
  else if ( distinct == 2 )
  {
    dif1 = y2;
    dif2 = ( ( y3 - y1 ) / ( x3 - x1 )
             - y2 ) / ( x3 - x2 );
  }
  else
  {
    dif1 = ( y2 - y1 ) / ( x2 - x1 );
    dif2 = ( ( y3 - y1 ) / ( x3 - x1 )
           - ( y2 - y1 ) / ( x2 - x1 ) ) / ( x3 - x2 );
  }



  *y = y1 + ( x - x1 ) * dif1 + ( x - x1 ) * ( x - x2 ) * dif2;
  *yp = dif1 + ( 2.0 * x - x1 - x2 ) * dif2;
  *ypp = 2.0 * dif2;

  return;
}


void r8poly2_val2 ( int ndata, double tdata[],
  double ydata[], int left, double tval, double *yval )

















































{
  double dif1;
  double dif2;
  double t1;
  double t2;
  double t3;
  double y1;
  double y2;
  double y3;



  if ( left < 0 || ndata-3 < left )
  {
    cerr << "\n";
    cerr << "RPOLY2_VAL2 - Fatal error!\n";
    cerr << "  LEFT < 0 or NDATA-3 < LEFT.\n";
    exit ( 1 );
  }



  t1 = tdata[left];
  t2 = tdata[left+1];
  t3 = tdata[left+2];

  if ( t2 <= t1 || t3 <= t2 )
  {
    cerr << "\n";
    cerr << "RPOLY2_VAL2 - Fatal error!\n";
    cerr << "  T2 <= T1 or T3 <= T2.\n";
    cerr << "  T1 = " << t1 << "\n";
    cerr << "  T2 = " << t2 << "\n";
    cerr << "  T3 = " << t3 << "\n";
    exit ( 1 );
  }



  y1 = ydata[left];
  y2 = ydata[left+1];
  y3 = ydata[left+2];

  dif1 = ( y2 - y1 ) / ( t2 - t1 );
  dif2 =
    ( ( y3 - y1 ) / ( t3 - t1 )
    - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 );

  *yval = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 );

  return;
}


void r8pp_delete ( int m, int n, double **a )




































{
  int i;

  for ( i = 0; i < m; i++ )
  {
    delete [] a[i];
  }

  delete [] a;

  return;
}


double **r8pp_new ( int m, int n )








































{
  double **a;
  int i;

  a = new double *[m];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "R8PP_NEW - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( i = 0; i < m; i++ )
  {
    a[i] = new double[n];
    if ( a[i] == NULL )
    {
      cerr << "\n";
      cerr << "R8PP_NEW - Fatal error!\n";
      cerr << "  Unable to allocate row array.\n";
      exit ( 1 );
    }
  }

  return a;
}



int r8r8_compare ( double x1, double y1, double x2, double y2 )


































{
  int value;

  if ( x1 < x2 )
  {
    value = -1;
  }
  else if ( x2 < x1 )
  {
    value = +1;
  }
  else if ( y1 < y2 )
  {
    value = -1;
  }
  else if ( y2 < y1 )
  {
    value = +1;
  }
  else
  {
    value = 0;
  }

  return value;
}


void r8r8_print ( double a1, double a2, string title )



































{
  cout << "  " << title << " : ";
  cout << "  ( " << setw(12) << a1
       << ", "   << setw(12) << a2 << " )\n";

  return;
}


int r8r8r8_compare ( double x1, double y1, double z1, double x2, double y2,
  double z2 )


































{
  int value;

  if ( x1 < x2 )
  {
    value = -1;
  }
  else if ( x2 < x1 )
  {
    value = +1;
  }
  else if ( y1 < y2 )
  {
    value = -1;
  }
  else if ( y2 < y1 )
  {
    value = +1;
  }
  else if ( z1 < z2 )
  {
    value = -1;
  }
  else if ( z2 < z1 )
  {
    value = +1;
  }
  else
  {
    value = 0;
  }

  return value;
}


void r8r8r8vec_index_insert_unique ( int maxn, int &n, double x[], double y[],
  double z[], int indx[], double xval, double yval, double zval, int &ival,
  int &ierror )










































{
  int equal;
  int i;
  int less;
  int more;

  ierror = 0;

  if ( n <= 0 )
  {
    if ( maxn <= 0 )
    {
      ierror = 1;
      cerr << "\n";
      cerr << "R8R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!\n";
      cerr << "  Not enough space to store new data.\n";
      return;
    }
    n = 1;
    x[0] = xval;
    y[0] = yval;
    z[0] = zval;
    indx[0] = 1;
    ival = 1;
    return;
  }



  r8r8r8vec_index_search ( n, x, y, z, indx, xval, yval, zval,
    less, equal, more );

  if ( equal == 0 )
  {
    if ( maxn <= n )
    {
      ierror = 1;
      cerr << "\n";
      cerr << "R8R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!\n";
      cerr << "  Not enough space to store new data.\n";
      return;
    }

    x[n] = xval;
    y[n] = yval;
    z[n] = zval;
    ival = n + 1;
    for ( i = n - 1; more - 1 <= i; i-- )
    {
      indx[i+1] = indx[i];
    }
    
    indx[more-1] = n + 1;
    n = n + 1;
  }
  else
  {
    ival = indx[equal-1];
  }

  return;
}


void r8r8r8vec_index_search ( int n, double x[], double y[], double z[],
  int indx[], double xval, double yval, double zval, int &less, int &equal,
  int &more )



































{
  int compare;
  int hi;
  int lo;
  int mid;
  double xhi;
  double xlo;
  double xmid;
  double yhi;
  double ylo;
  double ymid;
  double zhi;
  double zlo;
  double zmid;

  if ( n <= 0 )
  {
    less = 0;
    equal = 0;
    more = 0;
    return;
  }

  lo = 1;
  hi = n;

  xlo = x[indx[lo-1]-1];
  ylo = y[indx[lo-1]-1];
  zlo = z[indx[lo-1]-1];

  xhi = x[indx[hi-1]-1];
  yhi = y[indx[hi-1]-1];
  zhi = z[indx[hi-1]-1];

  compare = r8r8r8_compare ( xval, yval, zval, xlo, ylo, zlo );

  if ( compare == -1 )
  {
    less = 0;
    equal = 0;
    more = 1;
    return;
  }
  else if ( compare == 0 )
  {
    less = 0;
    equal = 1;
    more = 2;
    return;
  }

  compare = r8r8r8_compare ( xval, yval, zval, xhi, yhi, zhi );

  if ( compare == 1 )
  {
    less = n;
    equal = 0;
    more = n + 1;
    return;
  }
  else if ( compare == 0 )
  {
    less = n - 1;
    equal = n;
    more = n + 1;
    return;
  }

  for ( ; ; )
  {
    if ( lo + 1 == hi )
    {
      less = lo;
      equal = 0;
      more = hi;
      return;
    }

    mid = ( lo + hi ) / 2;
    xmid = x[indx[mid-1]-1];
    ymid = y[indx[mid-1]-1];
    zmid = z[indx[mid-1]-1];

    compare = r8r8r8_compare ( xval, yval, zval, xmid, ymid, zmid );

    if ( compare == 0 )
    {
      equal = mid;
      less = mid - 1;
      more = mid + 1;
      return;
    }
    else if ( compare == -1 )
    {
      hi = mid;
    }
    else if ( compare == +1 )
    {
      lo = mid;
    }
  }

  return;
}


void r8r8vec_index_insert_unique ( int maxn, int &n, double x[], double y[],
  int indx[], double xval, double yval, int &ival, int &ierror )










































{
  int equal;
  int i;
  int less;
  int more;

  ierror = 0;

  if ( n <= 0 )
  {
    if ( maxn <= 0 )
    {
      cerr << "\n";
      cerr << "R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!\n";
      cerr << "  Not enough space to store new data.\n";
      exit ( 1 );
    }

    n = 1;
    x[0] = xval;
    y[0] = yval;
    indx[0] = 1;
    ival = 1;
    return;
  }



  r8r8vec_index_search ( n, x, y, indx, xval, yval, less, equal, more );

  if ( equal == 0 )
  {
    if ( maxn <= n )
    {
      cerr << "\n";
      cerr << "R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!\n";
      cerr << "  Not enough space to store new data.\n";
      exit ( 1 );
    }

    x[n] = xval;
    y[n] = yval;
    ival = n + 1;
    for ( i = n - 1; more - 1 <= i; i-- )
    {
      indx[i+1] = indx[i];
    }
    indx[more-1] = n + 1;
    n = n + 1;
  }
  else
  {
    ival = indx[equal-1];
  }

  return;
}


void r8r8vec_index_search ( int n, double x[], double y[], int indx[],
  double xval, double yval, int &less, int &equal, int &more )



































{
  int compare;
  int hi;
  int lo;
  int mid;
  double xhi;
  double xlo;
  double xmid;
  double yhi;
  double ylo;
  double ymid;

  if ( n <= 0 )
  {
    less = 0;
    equal = 0;
    more = 0;
    return;
  }

  lo = 1;
  hi = n;

  xlo = x[indx[lo-1]-1];
  ylo = y[indx[lo-1]-1];

  xhi = x[indx[hi-1]-1];
  yhi = y[indx[hi-1]-1];

  compare = r8r8_compare ( xval, yval, xlo, ylo );

  if ( compare == -1 )
  {
    less = 0;
    equal = 0;
    more = 1;
    return;
  }
  else if ( compare == 0 )
  {
    less = 0;
    equal = 1;
    more = 2;
    return;
  }

  compare = r8r8_compare ( xval, yval, xhi, yhi );

  if ( compare == 1 )
  {
    less = n;
    equal = 0;
    more = n + 1;
    return;
  }
  else if ( compare == 0 )
  {
    less = n - 1;
    equal = n;
    more = n + 1;
    return;
  }

  for ( ; ; )
  {
    if ( lo + 1 == hi )
    {
      less = lo;
      equal = 0;
      more = hi;
      return;
    }

    mid = ( lo + hi ) / 2;
    xmid = x[indx[mid-1]-1];
    ymid = y[indx[mid-1]-1];

    compare = r8r8_compare ( xval, yval, xmid, ymid );

    if ( compare == 0 )
    {
      equal = mid;
      less = mid - 1;
      more = mid + 1;
      return;
    }
    else if ( compare == -1 )
    {
      hi = mid;
    }
    else if ( compare == +1 )
    {
      lo = mid;
    }
  }

  return;
}


double **r8rmat_copy_new ( int m, int n, double **a )










































{
  double **b;
  int i;
  int j;

  b = r8rmat_new ( m, n );

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      b[i][j] = a[i][j];
    }
  }
  return b;
}


void r8rmat_delete ( int m, int n, double **a )



































{
  int i;

  for ( i = 0; i < m; i++ )
  {
    delete [] a[i];
  }

  delete [] a;

  return;
}


double *r8rmat_fs_new ( int n, double **a, double b[] )






























{
  double **a2;
  int i;
  int j;
  int k;
  int p;
  double t;
  double *x;

  a2 = r8rmat_copy_new ( n, n, a );
  x = r8vec_copy_new ( n, b );

  for ( k = 0; k < n; k++ )
  {



    p = k;

    for ( i = k + 1; i < n; i++ )
    {
      if ( fabs ( a2[p][k] ) < fabs ( a2[i][k] ) )
      {
        p = i;
      }
    }

    if ( a2[p][k] == 0.0 )
    {
      cerr << "\n";
      cerr << "R8RMAT_FS_NEW - Fatal error!\n";
      cerr << "  Zero pivot on step " << k << "\n";
      exit ( 1 );
    }



    if ( k != p )
    {
      for ( j = 0; j < n; j++ )
      {
        t        = a2[k][j];
        a2[k][j] = a2[p][j];
        a2[p][j] = t;
      }
      t    = x[k];
      x[k] = x[p];
      x[p] = t;
    }



    t = a2[k][k];
    a2[k][k] = 1.0;
    for ( j = k + 1; j < n; j++ )
    {
      a2[k][j] = a2[k][j] / t;
    }
    x[k] = x[k] / t;



    for ( i = k + 1; i < n; i++ )
    {
      if ( a2[i][k] != 0.0 )
      {
        t = - a2[i][k];
        a2[i][k] = 0.0;
        for ( j = k + 1; j < n; j++ )
        {
          a2[i][j] = a2[i][j] + t * a2[k][j];
        }
        x[i] = x[i] + t * x[k];
      }
    }
  }



  for ( j = n - 1; 1 <= j; j-- )
  {
    for ( i = 0; i < j; i++ )
    {
      x[i] = x[i] - a2[i][j] * x[j];
    }
  }

  r8rmat_delete ( n, n, a2 );

  return x;
}


double **r8rmat_new ( int m, int n )







































{
  double **a;
  int i;

  a = new double *[m];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "R8RMAT_NEW - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( i = 0; i < m; i++ )
  {
    a[i] = new double[n];
    if ( a[i] == NULL )
    {
      cerr << "\n";
      cerr << "R8RMAT_NEW - Fatal error!\n";
      cerr << "  Unable to allocate row array.\n";
      exit ( 1 );
    }
  }

  return a;
}


void r8rmat_print ( int m, int n, double **a, string title )





































{
  r8rmat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}


void r8rmat_print_some ( int m, int n, double **a, int ilo, int jlo, int ihi,
  int jhi, string title )










































{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }



  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }
    cout << "\n";





    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j - 1 << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";



    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( ihi < m )
    {
      i2hi = ihi;
    }
    else
    {
      i2hi = m;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {



      cout << setw(5) << i - 1 << ": ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1][j-1] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}


double *r8rmat_to_r8mat ( int m, int n, double **a )







































{
  double *b;
  int i;
  int j;

  b = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      b[i+j*m] = a[i][j];
    }
  }

  return b;
}


double **r8rmat_zeros ( int m, int n )







































{
  double **a;
  int i;
  int j;

  a = new double *[m];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "R8RMAT_ZEROS - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( i = 0; i < m; i++ )
  {
    a[i] = new double[n];
    if ( a[i] == NULL )
    {
      cerr << "\n";
      cerr << "R8RMAT_ZEROS - Fatal error!\n";
      cerr << "  Unable to allocate row array.\n";
      exit ( 1 );
    }
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i][j] = 0.0;
    }
  }
  return a;
}


void r8slmat_print ( int m, int n, double a[], string title )








































{
  int i;
  int indx;
  int j;
  int jhi;
  int jlo;
  int jmax;
  int nn;

  cout << "\n";
  cout << title << "\n";

  jmax = i4_min ( n, m - 1 );

  nn = 5;

  for ( jlo = 1; jlo <= jmax; jlo = jlo + nn )
  {
    jhi = i4_min ( jlo + nn - 1, i4_min ( m - 1, jmax ) );
    cout << "\n";
    cout << "  Col   ";
    for ( j = jlo; j <= jhi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    for ( i = jlo + 1; i <= m; i++ )
    {
      cout << setw(5) << i << ":";
      jhi = i4_min ( jlo + nn - 1, i4_min ( i - 1, jmax ) );
      for ( j = jlo; j <= jhi; j++ )
      {
        indx = ( j - 1 ) * m + i - ( j * ( j + 1 ) ) / 2;
        cout << " " << setw(12) << a[indx-1];
      }
      cout << "\n";
    }
  }

  return;
}


void r8vec_01_to_ab ( int n, double a[], double amax, double amin )







































{
  double amax2;
  double amax3;
  double amin2;
  double amin3;
  int i;

  if ( amax == amin )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = amin;
    }
    return;
  }

  amax2 = r8_max ( amax, amin );
  amin2 = r8_min ( amax, amin );

  amin3 = r8vec_min ( n, a );
  amax3 = r8vec_max ( n, a );

  if ( amax3 != amin3 )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( amax3 - a[i]         ) * amin2
             + (         a[i] - amin3 ) * amax2 )
             / ( amax3          - amin3 );
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = 0.5 * ( amax2 + amin2 );
    }
  }

  return;
}


void r8vec_add ( int n, double a1[], double a2[] )
































{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a2[i] + a1[i];
  }
  return;
}


double r8vec_amax ( int n, double a[] )
































{
  double amax;
  int i;

  amax = 0.0;
  for ( i = 0; i < n; i++ )
  {
    if ( amax < fabs ( a[i] ) )
    {
      amax = fabs ( a[i] );
    }
  }

  return amax;
}


int r8vec_amax_index ( int n, double a[] )































{
  double amax;
  int amax_index;
  int i;

  if ( n <= 0 )
  {
    amax_index = -1;
  }
  else
  {
    amax_index = 1;
    amax = fabs ( a[0] );

    for ( i = 2; i <= n; i++ )
    {
      if ( amax < fabs ( a[i-1] ) )
      {
        amax_index = i;
        amax = fabs ( a[i-1] );
      }
    }
  }

  return amax_index;
}


double r8vec_amin ( int n, double a[] )
































{
  int i;
  const double r8_huge = 1.79769313486231571E+308;
  double value;

  value = r8_huge;
  for ( i = 0; i < n; i++ )
  {
    if ( fabs ( a[i] ) < value )
    {
      value = fabs ( a[i] );
    }
  }

  return value;
}


int r8vec_amin_index ( int n, double a[] )































{
  double amin;
  int amin_index;
  int i;

  if ( n <= 0 )
  {
    amin_index = -1;
  }
  else
  {
    amin_index = 1;
    amin = fabs ( a[0] );

    for ( i = 2; i <= n; i++ )
    {
      if ( fabs ( a[i-1] ) < amin )
      {
        amin_index = i;
        amin = fabs ( a[i-1] );
      }
    }
  }

  return amin_index;
}


double *r8vec_any_normal ( int dim_num, double v1[] )





































{
  int i;
  int j;
  int k;
  double *v2;
  double vj;
  double vk;

  if ( dim_num < 2 )
  {
    cerr << "\n";
    cerr << "R8VEC_ANY_NORMAL - Fatal error!\n";
    cerr << "  Called with DIM_NUM < 2.\n";
    exit ( 1 );
  }

  v2 = new double[dim_num];

  if ( r8vec_norm ( dim_num, v1 ) == 0.0 )
  {
    r8vec_zeros ( dim_num, v2 );
    v2[0] = 1.0;
    return v2;
  }







  j = -1;
  vj = 0.0;

  k = -1;
  vk = 0.0;

  for ( i = 0; i < dim_num; i++ )
  {
    if ( fabs ( vk ) < fabs ( v1[i] ) || k == -1 )
    {
      if ( fabs ( vj ) < fabs ( v1[i] ) || j == -1 )
      {
        k = j;
        vk = vj;
        j = i;
        vj = v1[i];
      }
      else
      {
        k = i;
        vk = v1[i];
      }
    }
  }




  r8vec_zeros ( dim_num, v2 );

  v2[j] = -vk / sqrt ( vk * vk + vj * vj );
  v2[k] =  vj / sqrt ( vk * vk + vj * vj );

  return v2;
}


void r8vec_append ( int *n, double **a, double value )





























{
  double *a_old;
  int i;



  a_old = *a;



  *a = new double[*n+1];



  for ( i = 0; i < *n; i++ )
  {
    (*a)[i] = a_old[i];
  }
  (*a)[*n] = value;



  *n = *n + 1;



  delete [] a_old;

  return;
}


double *r8vec_append_new ( int n, double a[], double value )





































{
  double *b;
  int i;

  b = new double[n+1];

  for ( i = 0; i < n; i++ )
  {
    b[i] = a[i];
  }
  b[n] = value;

  return b;
}


double r8vec_asum ( int n, double a[] )































{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + fabs ( a[i] );
  }
  return value;
}


void r8vec_bin ( int n, double x[], int bin_num, double bin_min, double bin_max,
  int bin[], double bin_limit[] )


























































{
  int i;
  int j;
  double t;

  if ( bin_max == bin_min )
  {
    cerr << "\n";
    cerr << "R8VEC_BIN - Fatal error!\n";
    cerr << "  BIN_MIN = BIN_MAX = " << bin_max << ".\n";
    exit ( 1 );
  }

  for ( i = 0; i <= bin_num + 1; i++ )
  {
    bin[i] = 0;
  }

  for ( i = 0; i < n; i++ )
  {
    t = ( x[i] - bin_min ) / ( bin_max - bin_min );

    if ( t < 0.0 )
    {
      j = 0;
    }
    else if ( 1.0 <= t )
    {
      j = bin_num + 1;
    }
    else
    {
      j = 1 + ( int ) ( ( double ) ( bin_num ) * t );
    }
    bin[j] = bin[j] + 1;
  }



  for ( i = 0; i <= bin_num; i++ )
  {
    bin_limit[i] = (   ( double ) ( bin_num - i ) * bin_min   
                     + ( double ) (           i ) * bin_max ) 
                     / ( double ) ( bin_num     );
  }

  return;
}


void r8vec_binary_next ( int n, double bvec[] )






















































{
  int i;

  for ( i = n - 1; 0 <= i; i-- )
  {
    if ( bvec[i] == 0.0 )
    {
      bvec[i] = 1.0;
      return;
    }
    bvec[i] = 0.0;
  }

  return;
}


void r8vec_bracket ( int n, double x[], double xval, int &left, int &right )


















































{
  int i;

  for ( i = 2; i <= n - 1; i++ )
  {
    if ( xval < x[i-1] )
    {
      left = i - 1;
      right = i;
      return;
    }

   }

  left = n - 1;
  right = n;

  return;
}


void r8vec_bracket2 ( int n, double x[], double xval, int start, int &left,
  int &right )
































































{
  int high;
  int low;



  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R8VEC_BRACKET2 - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( start < 1 || n < start )
  {
    start = ( n + 1 ) / 2;
  }



  if ( x[start-1] == xval )
  {
    left = start;
    right = start;
    return;
  }



  else if ( x[start-1] < xval )
  {



    if ( n < start + 1 )
    {
      left = start;
      right = -1;
      return;
    }



    else if ( xval == x[start] )
    {
      left = start + 1;
      right = start + 1;
      return;
    }



    else if ( xval < x[start] )
    {
      left = start;
      right = start + 1;
      return;
    }



    else if ( n < start + 2 )
    {
      left = start + 1;
      right = -1;
      return;
    }



    else if ( xval == x[start+1] )
    {
      left = start + 2;
      right = start + 2;
      return;
    }



    else if ( xval < x[start+1] )
    {
      left = start + 1;
      right = start + 2;
      return;
    }




    else
    {
      low = start + 2;
      high = n;

      r8vec_bracket ( high + 1 - low, x+low-1, xval, left, right );

      left = left + low - 1;
      right = right + low - 1;
    }
  }



  else if ( start == 1 )
  {
    left = -1;
    right = start;
    return;
  }



  else if ( xval == x[start-2] )
  {
    left = start - 1;
    right = start - 1;
    return;
  }



  else if ( x[start-2] <= xval )
  {
    left = start - 1;
    right = start;
    return;
  }




  else
  {
    low = 1;
    high = start - 1;
    r8vec_bracket ( high + 1 - low, x, xval, left, right );
  }

  return;
}


void r8vec_bracket3 ( int n, double t[], double tval, int &left )























































{
  int high;
  int low;
  int mid;



  if ( n < 2 )
  {
    cerr << "\n";
    cerr << "R8VEC_BRACKET3 - Fatal error!\n";
    cerr << "  N must be at least 2.\n";
    exit ( 1 );
  }



  if ( left < 0 || n - 2 < left )
  {
    left = ( n - 1 ) / 2;
  }




  if ( tval < t[left] )
  {
    if ( left == 0 )
    {
      return;
    }
    else if ( left == 1 )
    {
      left = 0;
      return;
    }
    else if ( t[left-1] <= tval )
    {
      left = left - 1;
      return;
    }
    else if ( tval <= t[1] )
    {
      left = 0;
      return;
    }



    low = 1;
    high = left - 2;

    for ( ; ; )
    {
      if ( low == high )
      {
        left = low;
        return;
      }

      mid = ( low + high + 1 ) / 2;

      if ( t[mid] <= tval )
      {
        low = mid;
      }
      else
      {
        high = mid - 1;
      }
    }
  }




  else if ( t[left+1] < tval )
  {
    if ( left == n - 2 )
    {
      return;
    }
    else if ( left == n - 3 )
    {
      left = left + 1;
      return;
    }
    else if ( tval <= t[left+2] )
    {
      left = left + 1;
      return;
    }
    else if ( t[n-2] <= tval )
    {
      left = n - 2;
      return;
    }



    low = left + 2;
    high = n - 3;

    for ( ; ; )
    {

      if ( low == high )
      {
        left = low;
        return;
      }

      mid = ( low + high + 1 ) / 2;

      if ( t[mid] <= tval )
      {
        low = mid;
      }
      else
      {
        high = mid - 1;
      }
    }
  }




  else
  {
  }

  return;
}


void r8vec_bracket4 ( int nt, double t[], int ns, double s[], int left[] )






















































{
  int high;
  int i;
  int low;
  int mid;



  if ( nt < 2 )
  {
    cerr << "\n";
    cerr << "R8VEC_BRACKET4 - Fatal error!\n";
    cerr << "  NT must be at least 2.\n";
    exit ( 1 );
  }

  for ( i = 0; i < ns; i++ )
  {
    left[i] = ( nt - 1 ) / 2;




    if ( s[i] < t[left[i]] )
    {
      if ( left[i] == 0 )
      {
        continue;
      }
      else if ( left[i] == 1 )
      {
        left[i] = 0;
        continue;
      }
      else if ( t[left[i]-1] <= s[i] )
      {
        left[i] = left[i] - 1;
        continue;
      }
      else if ( s[i] <= t[1] )
      {
        left[i] = 0;
        continue;
      }



      low = 1;
      high = left[i] - 2;

      for ( ; ; )
      {
        if ( low == high )
        {
          left[i] = low;
          break;
        }

        mid = ( low + high + 1 ) / 2;

        if ( t[mid] <= s[i] )
        {
          low = mid;
        }
        else
        {
          high = mid - 1;
        }
      }
    }




    else if ( t[left[i]+1] < s[i] )
    {
      if ( left[i] == nt - 2 )
      {
        continue;
      }
      else if ( left[i] == nt - 3 )
      {
        left[i] = left[i] + 1;
        continue;
      }
      else if ( s[i] <= t[left[i]+2] )
      {
        left[i] = left[i] + 1;
        continue;
      }
      else if ( t[nt-2] <= s[i] )
      {
        left[i] = nt - 2;
        continue;
      }



      low = left[i] + 2;
      high = nt - 3;

      for ( ; ; )
      {

        if ( low == high )
        {
          left[i] = low;
          break;
        }

        mid = ( low + high + 1 ) / 2;

        if ( t[mid] <= s[i] )
        {
          low = mid;
        }
        else
        {
          high = mid - 1;
        }
      }
    }



    else
    {
    }
  }
  return;
}


int r8vec_bracket5 ( int nd, double xd[], double xi )









































{
  int b;
  int l;
  int m;
  int r;

  if ( xi < xd[0] || xd[nd-1] < xi )
  {
    b = -1;
  }
  else
  {
    l = 0;
    r = nd - 1;

    while ( l + 1 < r )
    {
      m = ( l + r ) / 2;
      if ( xi < xd[m] )
      {
        r = m;
      }
      else
      {
        l = m;
      }
    }
    b = l;
  }

  return b;
}


int *r8vec_bracket6 ( int nd, double xd[], int ni, double xi[] )











































{
  int *b;
  int i;
  int l;
  int m;
  int r;

  b = new int[ni];

  for ( i = 0; i < ni; i++ )
  {
    if ( xi[i] < xd[0] || xd[nd-1] < xi[i] )
    {
      b[i] = -1;
    }
    else
    {
      l = 0;
      r = nd - 1;

      while ( l + 1 < r )
      {
        m = ( l + r ) / 2;
        if ( xi[i] < xd[m] )
        {
          r = m;
        }
        else
        {
          l = m;
        }
      }

      b[i] = l;
    }
  }

  return b;
}


double *r8vec_cheby_extreme_new ( int n, double a, double b )































{
  double c;
  int i;
  const double r8_pi = 3.141592653589793;
  double theta;
  double *x;

  x = new double[n];

  if ( n == 1 )
  {
    x[0] = ( a + b ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      theta = ( double ) ( n - i - 1 ) * r8_pi / ( double ) ( n - 1 );

      c = cos ( theta );

      if ( ( n % 2 ) == 1 )
      {
        if ( 2 * i + 1 == n )
        {
          c = 0.0;
        }
      }

      x[i] = ( ( 1.0 - c ) * a  
             + ( 1.0 + c ) * b ) 
             /   2.0;
    }
  }

  return x;
}


double *r8vec_cheby_zero_new ( int n, double a, double b )































{
  double c;
  int i;
  const double r8_pi = 3.141592653589793;
  double theta;
  double *x;

  x = new double[n];

  if ( n == 1 )
  {
    x[0] = ( a + b ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      theta = ( double ) ( 2 * ( n - i ) - 1 ) * r8_pi / ( double ) ( 2 * n );

      c = cos ( theta );

      if ( ( n % 2 ) == 1 )
      {
        if ( 2 * i + 1 == n )
        {
          c = 0.0;
        }
      }

      x[i] = ( ( 1.0 - c ) * a  
             + ( 1.0 + c ) * b ) 
             /   2.0;
    }
  }

  return x;
}


double *r8vec_cheby1space_new ( int n, double a, double b )
































{
  double c;
  int i;
  const double r8_pi = 3.141592653589793;
  double theta;
  double *x;

  x = new double[n];

  if ( n == 1 )
  {
    x[0] = ( a + b ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      theta = ( double ) ( n - i - 1 ) * r8_pi / ( double ) ( n - 1 );

      c = cos ( theta );

      if ( ( n % 2 ) == 1 )
      {
        if ( 2 * i + 1 == n )
        {
          c = 0.0;
        }
      }

      x[i] = ( ( 1.0 - c ) * a  
             + ( 1.0 + c ) * b ) 
             /   2.0;
    }
  }

  return x;
}


double *r8vec_cheby2space_new ( int n, double a, double b )
































{
  double c;
  int i;
  const double r8_pi = 3.141592653589793;
  double theta;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    theta = ( double ) ( n - i ) * r8_pi / ( double ) ( n + 1 );

    c = cos ( theta );

    x[i] = ( ( 1.0 - c ) * a  
           + ( 1.0 + c ) * b ) 
           /   2.0;
  }

  return x;
}


int r8vec_compare ( int n, double a[], double b[] )















































{
  int isgn;
  int k;

  isgn = 0;

  for ( k = 0; k < n; k++ )
  {
    if ( a[k] < b[k] )
    {
      isgn = -1;
      return isgn;
    }
    else if ( b[k] < a[k] )
    {
      isgn = +1;
      return isgn;
    }
  }
  return isgn;
}


void r8vec_concatenate ( int n1, double a[], int n2, double b[], double c[] )



































{
  int i;

  for ( i = 0; i < n1; i++ )
  {
    c[i] = a[i];
  }
  for ( i = 0; i < n2; i++ )
  {
    c[n1+i] = b[i];
  }

  return;
}


double *r8vec_concatenate_new ( int n1, double a[], int n2, double b[] )



































{
  int i;
  double *c;

  c = new double[n1+n2];

  for ( i = 0; i < n1; i++ )
  {
    c[i] = a[i];
  }
  for ( i = 0; i < n2; i++ )
  {
    c[n1+i] = b[i];
  }

  return c;
}


double *r8vec_convolution ( int m, double x[], int n, double y[] )

































































{
  int i;
  int j;
  double *z;

  z = new double[m+n-1];

  for ( i = 0; i < m + n - 1; i++ )
  {
    z[i] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      z[j+i] = z[j+i] + x[i] * y[j];
    }
  }
  return z;
}


double *r8vec_convolution_circ ( int n, double x[], double y[] )


































































{
  int i;
  int m;
  double *z;

  z = new double[n];

  for ( m = 1; m <= n; m++ )
  {
    z[m-1] = 0.0;
    for ( i = 1; i <= m; i++ )
    {
      z[m-1] = z[m-1] + x[i-1] * y[m-i];
    }
    for ( i = m+1; i <= n; i++ )
    {
      z[m-1] = z[m-1] + x[i-1] * y[n+m-i];
    }
  }

  return z;
}


void r8vec_copy ( int n, double a1[], double a2[] )































{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}


double *r8vec_copy_new ( int n, double a1[] )































{
  double *a2;
  int i;

  a2 = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
}


double r8vec_correlation ( int n, double x[], double y[] )





































{
  double correlation;
  double x_norm;
  double xy_dot;
  double y_norm;

  x_norm = r8vec_norm ( n, x );
  y_norm = r8vec_norm ( n, y );
  xy_dot = r8vec_dot_product ( n, x, y );

  if ( x_norm == 0.0 || y_norm == 0.0 )
  {
    correlation = 0.0;
  }
  else
  {
    correlation = xy_dot / x_norm / y_norm;
  }

  return correlation;
}


double r8vec_covar ( int n, double x[], double y[] )



























{
  int i;
  double value;
  double x_average;
  double y_average;

  x_average = 0.0;
  for ( i = 0; i < n; i++ )
  {
    x_average = x_average + x[i];
  }
  x_average = x_average / ( double ) ( n );

  y_average = 0.0;
  for ( i = 0; i < n; i++ )
  {
    y_average = y_average + x[i];
  }
  y_average = y_average / ( double ) ( n );

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + ( x[i] - x_average ) * ( y[i] - y_average );
  }

  value = value / ( double ) ( n - 1 );

  return value;
}


double r8vec_cross_product_2d ( double v1[2], double v2[2] )































{
  double value;

  value = v1[0] * v2[1] - v1[1] * v2[0];

  return value;
}


double r8vec_cross_product_affine_2d ( double v0[2], double v1[2],
  double v2[2] )

































{
  double value;

  value =
      ( v1[0] - v0[0] ) * ( v2[1] - v0[1] )
    - ( v2[0] - v0[0] ) * ( v1[1] - v0[1] );

  return value;
}


double *r8vec_cross_product_3d ( double v1[3], double v2[3] )





























{
  double *v3;

  v3 = new double[3];

  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return v3;
}


double *r8vec_cross_product_affine_3d ( double v0[3], double v1[3],
  double v2[3] )































{
  double *v3;

  v3 = ( double * ) malloc ( 3 * sizeof ( double ) );

  v3[0] =
      ( v1[1] - v0[1] ) * ( v2[2] - v0[2] )
    - ( v2[1] - v0[1] ) * ( v1[2] - v0[2] );

  v3[1] =
      ( v1[2] - v0[2] ) * ( v2[0] - v0[0] )
    - ( v2[2] - v0[2] ) * ( v1[0] - v0[0] );

  v3[2] =
      ( v1[0] - v0[0] ) * ( v2[1] - v0[1] )
    - ( v2[0] - v0[0] ) * ( v1[1] - v0[1] );

  return v3;
}


double *r8vec_cum_new ( int n, double a[] )







































{
  double *a_cum;
  int i;

  a_cum = new double[n];

  a_cum[0] = a[0];

  for ( i = 1; i < n; i++ )
  {
    a_cum[i] = a_cum[i-1] + a[i];
  }

  return a_cum;
}


double *r8vec_cum0_new ( int n, double a[] )







































{
  double *a_cum;
  int i;

  a_cum = new double[n+1];

  a_cum[0] = 0.0;

  for ( i = 1; i <= n; i++ )
  {
    a_cum[i] = a_cum[i-1] + a[i-1];
  }

  return a_cum;
}


double *r8vec_dif ( int n, double h )





















































































{
  double *cof;
  int i;
  int j;

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_DIF - Fatal error!\n";
    cerr << "  Derivative order N = " << n << "\n";
    cerr << "  but N must be at least 0.\n";
    exit ( 1 );
  }

  if ( h <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8VEC_DIF - Fatal error!\n";
    cerr << "  The half sampling spacing is H = " << h << "\n";
    cerr << "  but H must be positive.\n";
    exit ( 1 );
  }

  cof = new double[n+1];

  for ( i = 0; i <= n; i++ )
  {
    cof[i] = 1.0;

    for ( j = i - 1; 1 <= j; j-- )
    {
      cof[j] = -cof[j] + cof[j-1];
    }

    if ( 0 < i )
    {
      cof[0] = - cof[0];
    }
  }

  for ( i = 0; i <= n; i++ )
  {
    cof[i] = cof[i] / pow ( 2.0 * h, n );
  }

  return cof;
}


double r8vec_diff_norm ( int n, double a[], double b[] )



































{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( a[i] - b[i] ) * ( a[i] - b[i] );
  }
  value = sqrt ( value );

  return value;
}


double r8vec_diff_norm_l1 ( int n, double a[], double b[] )



































{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + fabs ( a[i] - b[i] );
  }
  return value;
}


double r8vec_diff_norm_l2 ( int n, double a[], double b[] )



































{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( a[i] - b[i] ) * ( a[i] - b[i] );
  }
  value = sqrt ( value );

  return value;
}


double r8vec_diff_norm_li ( int n, double a[], double b[] )



































{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = r8_max ( value, fabs ( a[i] - b[i] ) );
  }
  return value;
}


double r8vec_diff_norm_squared ( int n, double a[], double b[] )



































{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( a[i] - b[i] ) * ( a[i] - b[i] );
  }

  return value;
}


void r8vec_direct_product ( int factor_index, int factor_order,
  double factor_value[], int factor_num, int point_num, double x[] )

























































































































{
  static int contig = 0;
  int i;
  int j;
  int k;
  static int rep = 0;
  static int skip = 0;
  int start;

  if ( factor_index == 0 )
  {
    contig = 1;
    skip = 1;
    rep = point_num;
    for ( j = 0; j < point_num; j++ )
    {
      for ( i = 0; i < factor_num; i++ )
      {
        x[i+j*factor_num] = 0.0;
      }
    }
  }

  rep = rep / factor_order;
  skip = skip * factor_order;

  for ( i = 0; i < factor_order; i++ )
  {
    start = 0 + i * contig;

    for ( k = 1; k <= rep; k++ )
    {
      for ( j = start; j < start + contig; j++ )
      {
        x[factor_index+j*factor_num] = factor_value[i];
      }
      start = start + skip;
    }
  }
  contig = contig * factor_order;

  return;
}


void r8vec_direct_product2 ( int factor_index, int factor_order,
  double factor_value[], int factor_num, int point_num, double w[] )

























































































































{
  static int contig = 0;
  int i;
  int j;
  int k;
  static int rep = 0;
  static int skip = 0;
  int start;

  if ( factor_index == 0 )
  {
    contig = 1;
    skip = 1;
    rep = point_num;
    for ( i = 0; i < point_num; i++ )
    {
      w[i] = 1.0;
    }
  }

  rep = rep / factor_order;
  skip = skip * factor_order;

  for ( j = 0; j < factor_order; j++ )
  {
    start = 0 + j * contig;

    for ( k = 1; k <= rep; k++ )
    {
      for ( i = start; i < start + contig; i++ )
      {
        w[i] = w[i] * factor_value[j];
      }
      start = start + skip;
    }
  }

  contig = contig * factor_order;

  return;
}


double r8vec_distance ( int dim_num, double v1[], double v2[] )
































{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < dim_num; i++ )
  {
    value = pow ( v1[i] - v2[i], 2 );
  }
  value = sqrt ( value );

  return value;
}


void r8vec_divide ( int n, double a[], double s )
































{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] / s;
  }
  return;
}


double r8vec_dot_product ( int n, double a1[], double a2[] )































{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}


double r8vec_dot_product_affine ( int n, double v0[], double v1[], double v2[] )





























{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + ( v1[i] - v0[i] ) * ( v2[i] - v0[i] );
  }
  return value;
}


double r8vec_entropy ( int n, double x[] )




































{
  int i;
  double value;
  double x_sum;
  double xi;

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] < 0.0 )
    {
      cerr << "\n";
      cerr << "R8VEC_ENTROPY - Fatal error!\n";
      cerr << "  Some entries are negative.\n";
      exit ( 1 );
    }
  }

  x_sum = 0.0;
  for ( i = 0; i < n; i++ )
  {
    x_sum = x_sum + x[i];
  }

  if ( x_sum == 0.0 )
  {
    cerr << "\n";
    cerr << "R8VEC_ENTROPY - Fatal error!\n";
    cerr << "  Entries sum to 0.\n";
    exit ( 1 );
  }

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    if ( 0.0 < x[i] )
    {
      xi = x[i] / x_sum;
      value = value - r8_log_2 ( xi ) * xi;
    }
  }

  return value;
}


bool r8vec_eq ( int n, double a1[], double a2[] )
































{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a1[i] != a2[i] )
    {
      return false;
    }
  }
  return true;
}


void r8vec_even ( int n, double alo, double ahi, double a[] )

































{
  int i;

  if ( n == 1 )
  {
    a[0] = 0.5 * ( alo + ahi );
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - i - 1 ) * alo
             + ( double ) (     i     ) * ahi )
             / ( double ) ( n     - 1 );
    }
  }

  return;
}


double *r8vec_even_new ( int n, double alo, double ahi )

































{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = 0.5 * ( alo + ahi );
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - i - 1 ) * alo
             + ( double ) (     i     ) * ahi )
             / ( double ) ( n     - 1 );
    }
  }

  return a;
}


double r8vec_even_select ( int n, double xlo, double xhi, int ival )







































{
  double xval;

  if ( n == 1 )
  {
    xval = 0.5 * ( xlo + xhi );
  }
  else
  {
    xval = ( ( double ) ( n - ival     ) * xlo
           + ( double ) (     ival - 1 ) * xhi )
           / ( double ) ( n        - 1 );
  }

  return xval;
}


void r8vec_even2 ( int maxval, int nfill[], int nold, double xold[],
  int &nval, double xval[] )


























































{
  int i;
  int j;
  int nadd;

  nval = 1;

  for ( i = 1; i <= nold - 1; i++ )
  {

    if ( nfill[i-1] < 0 )
    {
      cerr << "\n";
      cerr << "R8VEC_EVEN2 - Fatal error!\n";
      cerr << "  NFILL[I-1] is negative for I = " << i << "\n";
      cerr << "  NFILL[I-1] = " << nfill[i-1] << "\n";
      exit ( 1 );
    }

    if ( maxval < nval + nfill[i-1] + 1 )
    {
      cerr << "\n";
      cerr << "R8VEC_EVEN2 - Fatal error!\n";
      cerr << "  MAXVAL = " << maxval << " is not large enough.\n";
      cerr << "  for the storage for interval I = " << i << "\n";
      exit ( 1 );
    }

    nadd = nfill[i-1] + 2;

    for ( j = 1; j <= nadd; j++ )
    {
      xval[nval+j-2] = ( ( double ) ( nadd - j     ) * xold[i-1]
                       + ( double ) (        j - 1 ) * xold[i] )
                       / ( double ) ( nadd     - 1 );
    }

    nval = nval + nfill[i-1] + 1;
  }

  return;
}


double r8vec_even2_select ( int n, double xlo, double xhi, int ival )










































{
  double xval;

  xval = ( ( double ) ( 2 * n - 2 * ival + 1 ) * xlo
         + ( double ) (         2 * ival - 1 ) * xhi )
         / ( double ) ( 2 * n                );

  return xval;
}


void r8vec_even3 ( int nold, int nval, double xold[], double xval[] )

















































{
  double density;
  int i;
  int ival;
  int j;
  int nmaybe;
  int npts;
  int ntemp;
  int ntot;
  double xlen;
  double xleni;
  double xlentot;

  xlen = 0.0;
  for ( i = 1; i <= nold - 1; i++ )
  {
    xlen = xlen + fabs ( xold[i] - xold[i-1] );
  }

  ntemp = nval - nold;

  density = ( double ) ( ntemp ) / xlen;

  ival = 1;
  ntot = 0;
  xlentot = 0.0;

  for ( i = 1; i <= nold - 1; i++ )
  {
    xleni = fabs ( xold[i] - xold[i-1] );
    npts = ( int ) ( density * xleni );
    ntot = ntot + npts;





    xlentot = xlentot + xleni;
    nmaybe = r8_nint ( xlentot * density );

    if ( ntot < nmaybe )
    {
      npts = npts + nmaybe - ntot;
      ntot = nmaybe;
    }
    for ( j = 1; j <= npts + 2; j++ )
    {
      xval[ival+j-2] = ( ( double ) ( npts+2 - j     ) * xold[i-1]
                       + ( double ) (          j - 1 ) * xold[i] )
                       / ( double ) ( npts+2     - 1 );
    }
    ival = ival + npts + 1;
  }

  return;
}


double *r8vec_expand_linear ( int n, double x[], int fat )


































{
  int i;
  int j;
  int k;
  double *xfat;

  xfat = new double[(n-1)*(fat+1)+1];

  k = 0;

  for ( i = 0; i < n-1; i++ )
  {
    xfat[k] = x[i];
    k = k + 1;

    for ( j = 1; j <= fat; j++ )
    {
      xfat[k] = ( ( double ) ( fat - j + 1 ) * x[i]
                + ( double ) (       j     ) * x[i+1] )
                / ( double ) ( fat     + 1 );
      k = k + 1;
    }
  }

  xfat[k] = x[n-1];
  k = k + 1;

  return xfat;
}


double *r8vec_expand_linear2 ( int n, double x[], int before, int fat, 
  int after )
































































{
  int i;
  int j;
  int k;
  double *xfat;

  xfat = new double[before+(n-1)*(fat+1)+1+after];

  k = 0;



  for ( j = 1 - before + fat; j <= fat; j++ )
  {
    xfat[k] = ( ( double ) ( fat - j + 1 ) * ( x[0] - ( x[1] - x[0] ) ) 
              + ( double ) (       j     ) *   x[0]          ) 
              / ( double ) ( fat     + 1 );
    k = k + 1;
  }



  for ( i = 0; i < n - 1; i++ )
  {
    xfat[k] = x[0];
    k = k + 1;
    for ( j = 1; j <= fat; j++ )
    {
      xfat[k] = ( ( double ) ( fat - j + 1 ) * x[i]
                + ( double ) (       j     ) * x[i+1] ) 
                / ( double ) ( fat     + 1 );
      k = k + 1;
    }
  }

  xfat[k] = x[n-1];
  k = k + 1;



  for ( j = 1; j <= after; j++ )
  {
    xfat[k] = ( ( double ) ( fat - j + 1 ) * x[n-1]
              + ( double ) (       j     ) * ( x[n-1] + ( x[n-1] - x[n-2] ) ) ) 
              / ( double ) ( fat     + 1 );
    k = k + 1;
  }

  return xfat;
}


void r8vec_fill ( int n, double value, double x[] )



























{
  int i;

  for ( i = 0; i < n; i++ )
  {
    x[i] = value;
  }
  return;
}


double *r8vec_fill_new ( int n, double value )



























{
  int i;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = value;
  }
  return x;
}


int *r8vec_first_index ( int n, double a[], double tol )




































{
  int *first_index;
  int i;
  int j;

  first_index = new int[n];

  for ( i = 0; i < n; i++ )
  {
    first_index[i] = -1;
  }
  for ( i = 0; i < n; i++ )
  {
    if ( first_index[i] == -1 )
    {
      first_index[i] = i;
      for ( j = i + 1; j < n; j++ )
      {
        if ( fabs ( a[i] - a[j] ) <= tol )
        {
          first_index[j] = i;
        }
      }
    }
  }
  return first_index;
}


double r8vec_frac ( int n, double a[], int k )




































{
  double frac;
  int i;
  int iryt;
  int j;
  int left;
  double temp;
  double x;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_FRAC - Fatal error!\n";
    cerr << "  Illegal nonpositive value of N = " << n << "\n";
    exit ( 1 );
  }

  if ( k <= 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_FRAC - Fatal error!\n";
    cerr << "  Illegal nonpositive value of K = " << k << "\n";
    exit ( 1 );
  }

  if ( n < k )
  {
    cerr << "\n";
    cerr << "R8VEC_FRAC - Fatal error!\n";
    cerr << "  Illegal N < K, K = " << k << "\n";
    exit ( 1 );
  }

  left = 1;
  iryt = n;

  for ( ; ; )
  {
    if ( iryt <= left )
    {
      frac = a[k-1];
      break;
    }

    x = a[k-1];
    i = left;
    j = iryt;

    for ( ; ; )
    {
      if ( j < i )
      {
        if ( j < k )
        {
          left = i;
        }
        if ( k < i )
        {
          iryt = j;
        }
        break;
      }



      while ( a[i-1] < x )
      {
        i = i + 1;
      }



      while ( x < a[j-1] )
      {
        j = j - 1;
      }

      if ( i <= j )
      {
        temp   = a[i-1];
        a[i-1] = a[j-1];
        a[j-1] = temp;
        i = i + 1;
        j = j - 1;
      }
    }
  }

  return frac;
}


double *r8vec_fraction ( int n, double x[] )





















































{
  double *fraction;
  int i;

  fraction = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fraction[i] = fabs ( x[i] ) - ( double ) ( ( int ) ( fabs ( x[i] ) ) );
  }

  return fraction;
}


bool r8vec_gt ( int n, double a1[], double a2[] )






































{
  int i;

  for ( i = 0; i < n; i++ )
  {

    if ( a2[i] < a1[i] )
    {
       return true;
    }
    else if ( a1[i] < a2[i] )
    {
      return false;
    }

  }

  return false;
}


void r8vec_heap_a ( int n, double a[] )




















































{
  int i;
  int ifree;
  double key;
  int m;



  for ( i = (n/2)-1; 0 <= i; i-- )
  {




    key = a[i];
    ifree = i;

    for ( ; ; )
    {




      m = 2 * ifree + 1;



      if ( n <= m )
      {
        break;
      }
      else
      {



        if ( m + 1 < n )
        {




          if ( a[m+1] < a[m] )
          {
            m = m + 1;
          }
        }





        if ( a[m] <= key )
        {
          break;
        }
        a[ifree] = a[m];
        ifree = m;
      }
    }




    a[ifree] = key;
  }

  return;
}


void r8vec_heap_d ( int n, double a[] )




















































{
  int i;
  int ifree;
  double key;
  int m;



  for ( i = (n/2)-1; 0 <= i; i-- )
  {




    key = a[i];
    ifree = i;

    for ( ; ; )
    {




      m = 2 * ifree + 1;



      if ( n <= m )
      {
        break;
      }
      else
      {



        if ( m + 1 < n )
        {




          if ( a[m] < a[m+1] )
          {
            m = m + 1;
          }
        }





        if ( key < a[m] )
        {
          a[ifree] = a[m];
          ifree = m;
        }
        else
        {
          break;
        }
      }
    }




    a[ifree] = key;
  }

  return;
}


int *r8vec_histogram ( int n, double a[], double a_lo, double a_hi,
  int histo_num )









































{
  double delta;
  int *histo_gram;
  int i;
  int j;

  histo_gram = new int[histo_num+2];

  i4vec_zeros ( histo_num+2, histo_gram );

  delta = ( a_hi - a_lo ) / ( double ) ( 2 * histo_num );

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < a_lo )
    {
      histo_gram[0] = histo_gram[0] + 1;
    }
    else if ( a[i] <= a_hi )
    {
      j = r8_nint (
        ( ( a_hi -       delta - a[i]        ) * ( double ) ( 1         )
        + (      -       delta + a[i] - a_lo ) * ( double ) ( histo_num ) )
        / ( a_hi - 2.0 * delta        - a_lo ) );

      histo_gram[j] = histo_gram[j] + 1;
    }
    else if ( a_hi < a[i] )
    {
      histo_gram[histo_num+1] = histo_gram[histo_num+1] + 1;
    }
  }

  return histo_gram;
}


double *r8vec_house_column ( int n, double a_vec[], int k )









































{
  int i;
  double s;
  double *v;

  v = r8vec_zeros_new ( n );

  if ( k < 1 || n <= k )
  {
    return v;
  }

  s = r8vec_norm_l2 ( n+1-k, a_vec+k-1 );

  if ( s == 0.0 )
  {
    return v;
  }

  v[k-1] = a_vec[k-1] + fabs ( s ) * r8_sign ( a_vec[k-1] );

  r8vec_copy ( n-k, a_vec+k, v+k );



  s = r8vec_norm_l2 ( n-k+1, v+k-1 );

  for ( i = k - 1; i < n; i++ )
  {
    v[i] = v[i] / s;
  }

  return v;
}


double r8vec_i4vec_dot_product ( int n, double r8vec[], int i4vec[] )



































{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + r8vec[i] * ( double ) ( i4vec[i] );
  }
  return value;
}


double *r8vec_identity_row_new ( int n, int i )



























{
  double *a;
  int j;

  a = new double[n];

  for ( j = 0; j < n; j++ )
  {
    a[j] = 0.0;
  }

  if ( 0 <= i && i < n )
  {
    a[i] = 1.0;
  }

  return a;
}


void r8vec_index_delete_all ( int n, double x[], int indx[], double xval,
  int &n2, double x2[], int indx2[] )









































{
  int equal;
  int equal1;
  int equal2;
  int get;
  int i;
  int less;
  int more;
  int put;

  if ( n < 1 )
  {
    n2 = 0;
    return;
  }

  i4vec_copy ( n, indx, indx2 );
  r8vec_copy ( n, x, x2 );
  n2 = n;

  r8vec_index_search ( n2, x2, indx2, xval, less, equal, more );

  if ( equal == 0 )
  {
    return;
  }

  equal1 = equal;

  for ( ; ; )
  {
    if ( equal1 <= 1 )
    {
      break;
    }

    if ( x2[indx2[equal1-2]-1] != xval )
    {
      break;
    }
    equal1 = equal1 - 1;
  }

  equal2 = equal;

  for ( ; ; )
  {
    if ( n2 <= equal2 )
    {
      break;
    }

    if ( x2[indx2[equal2]-1] != xval )
    {
      break;
    }
    equal2 = equal2 + 1;
  }



  put = 0;

  for ( get = 1; get <= n2; get++ )
  {
    if ( x2[get-1] != xval )
    {
      put = put + 1;
      x2[put-1] = x2[get-1];
    }
  }



  for ( equal = equal1; equal <= equal2; equal++ )
  {
    for ( i = 1; i <= n2; i++ )
    {
      if ( indx2[equal-1] < indx2[i-1] )
      {
        indx2[i-1] = indx2[i-1] - 1;
      }
    }
  }



  for ( i = 0; i <= n2 - equal2 - 1; i++ )
  {
    indx2[equal1+i-1] = indx2[equal2+i];
  }
  for ( i = n2 + equal1 - equal2; i <= n2; i++ )
  {
    indx2[i-1] = 0;
  }



  n2 = put;

  return;
}


void r8vec_index_delete_dupes ( int n, double x[], int indx[],
  int &n2, double x2[], int indx2[] )














































{
  int i;
  int n3;
  double *x3;

  i = 0;
  n3 = 0;
  x3 = new double[n];

  for ( ; ; )
  {
    i = i + 1;

    if ( n < i )
    {
      break;
    }

    if ( 1 < i )
    {
      if ( x[indx[i-1]-1] == x3[n3-1] )
      {
        continue;
      }
    }
    n3 = n3 + 1;
    x3[n3-1] = x[indx[i-1]-1];
  }



  n2 = n3;
  r8vec_copy ( n3, x3, x2 );
  for ( i = 0; i < n3; i++ )
  {
    indx2[i] = i + 1;
  }

  delete [] x3;

  return;
}


void r8vec_index_delete_one ( int n, double x[], int indx[], double xval,
  int &n2, double x2[], int indx2[] )











































{
  int equal;
  int i;
  int j;
  int less;
  int more;

  if ( n < 1 )
  {
    n2 = 0;
    return;
  }

  n2 = n;
  i4vec_copy ( n2, indx, indx2 );
  r8vec_copy ( n2, x, x2 );

  r8vec_index_search ( n2, x2, indx2, xval, less, equal, more );

  if ( equal != 0 )
  {
    j = indx2[equal-1];
    for ( i = j; i <= n2 - 1; i++ )
    {
      x2[i-1] = x[i];
    }
    for ( i = equal; i <= n2 - 1; i++ )
    {
      indx2[i-1] = indx2[i];
    }
    for ( i = 1; i <= n2 - 1; i++ )
    {
      if ( j < indx2[i-1] )
      {
        indx2[i-1] = indx2[i-1] - 1;
      }
    }
    n2 = n2 - 1;
  }

  return;
}


void r8vec_index_insert ( int &n, double x[], int indx[], double xval )

































{
  int equal;
  int i;
  int less;
  int more;

  if ( n <= 0 )
  {
    n = 1;
    x[0] = xval;
    indx[0] = 1;
    return;
  }

  r8vec_index_search ( n, x, indx, xval, less, equal, more );

  x[n] = xval;
  for ( i = n; more <= i; i-- )
  {
    indx[i] = indx[i-1];
  }
  indx[more-1] = n + 1;
  n = n + 1;

  return;
}


void r8vec_index_insert_unique ( int &n, double x[], int indx[], double xval )






































{
  int equal;
  int i;
  int less;
  int more;

  if ( n <= 0 )
  {
    n = 1;
    x[0] = xval;
    indx[0] = 1;
    return;
  }



  r8vec_index_search ( n, x, indx, xval, less, equal, more );

  if ( equal == 0 )
  {
    x[n] = xval;
    for ( i = n; more <= i; i-- )
    {
      indx[i] = indx[i-1];
    }
    indx[more-1] = n + 1;
    n = n + 1;
  }

  return;
}


void r8vec_index_order ( int n, double x[], int indx[] )



































{
  int i;
  double *y;

  y = new double[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = x[indx[i]-1];
  }
  for ( i = 0; i < n; i++ )
  {
    x[i] = y[i];
  }
  delete [] y;

  return;
}


void r8vec_index_search ( int n, double x[], int indx[], double xval, int &less,
  int &equal, int &more )







































{
  int hi;
  int lo;
  int mid;
  double xhi;
  double xlo;
  double xmid;

  if ( n <= 0 )
  {
    less = 0;
    equal = 0;
    more = 0;
    return;
  }

  lo = 1;
  hi = n;
  xlo = x[indx[lo-1]-1];
  xhi = x[indx[hi-1]-1];

  if ( xval < xlo )
  {
    less = 0;
    equal = 0;
    more = 1;
    return;
  }
  else if ( xval == xlo )
  {
    less = 0;
    equal = 1;
    more = 2;
    return;
  }

  if ( xhi < xval )
  {
    less = n;
    equal = 0;
    more = n + 1;
    return;
  }
  else if ( xval == xhi )
  {
    less = n - 1;
    equal = n;
    more = n + 1;
    return;
  }

  for ( ; ; )
  {
    if ( lo + 1 == hi )
    {
      less = lo;
      equal = 0;
      more = hi;
      return;
    }

    mid = ( lo + hi ) / 2;
    xmid = x[indx[mid-1]-1];

    if ( xval == xmid )
    {
      equal = mid;
      less = mid - 1;
      more = mid + 1;
      return;
    }
    else if ( xval < xmid )
    {
      hi = mid;
    }
    else if ( xmid < xval )
    {
      lo = mid;
    }
  }
  return;
}


void r8vec_index_sort_unique ( int n, double x[], int &n2, double x2[],
  int indx2[] )



































{
  int i;

  n2 = 0;

  for ( i = 0; i < n; i++ )
  {
    r8vec_index_insert_unique ( n2, x2, indx2, x[i] );
  }

  for ( i = n2; i < n; i++ )
  {
    x2[i] = -1;
  }
  for ( i = n2; i < n; i++ )
  {
    indx2[i] = -1;
  }

  return;
}


void r8vec_index_sorted_range ( int n, double r[], int indx[], double r_lo,
  double r_hi, int &i_lo, int &i_hi )


































{
  int i1;
  int i2;
  int j1;
  int j2;



  if ( r[indx[n-1]] < r_lo )
  {
    i_lo = n;
    i_hi = n - 1;
    return;
  }

  if ( r_hi < r[indx[0]] )
  {
    i_lo = 0;
    i_hi = -1;
    return;
  }



  if ( n == 1 )
  {
    if ( r_lo <= r[indx[0]] && r[indx[0]] <= r_hi )
    {
      i_lo = 0;
      i_hi = 0;
    }
    else
    {
      i_lo = -1;
      i_hi = -2;
    }
    return;
  }



  if ( r_lo <= r[indx[0]] )
  {
    i_lo = 0;
  }
  else
  {





    j1 = 0;
    j2 = n - 1;
    i1 = ( j1 + j2 - 1 ) / 2;
    i2 = i1 + 1;

    for ( ; ; )
    {
      if ( r_lo < r[indx[i1]] )
      {
        j2 = i1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else if ( r[indx[i2]] < r_lo )
      {
        j1 = i2;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else
      {
        i_lo = i1;
        break;
      }
    }
  }



  if ( r[indx[n-1]] <= r_hi )
  {
    i_hi = n - 1;
  }
  else
  {
    j1 = i_lo;
    j2 = n - 1;
    i1 = ( j1 + j2 - 1 ) / 2;
    i2 = i1 + 1;

    for ( ; ; )
    {
      if ( r_hi < r[indx[i1]] )
      {
        j2 = i1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else if ( r[indx[i2]] < r_hi )
      {
        j1 = i2;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else
      {
        i_hi = i2;
        break;
      }
    }
  }







  if ( r[indx[i_lo]] < r_lo )
  {
    i_lo = i_lo + 1;
    if ( n - 1 < i_lo )
    {
      i_hi = i_lo - 1;
    }
  }

  if ( r_hi < r[indx[i_hi]] )
  {
    i_hi = i_hi - 1;
    if ( i_hi < 0 )
    {
      i_lo = i_hi + 1;
    }
  }

  return;
}


void r8vec_indexed_heap_d ( int n, double a[], int indx[] )


















































{
  int i;
  int ifree;
  int key;
  int m;



  for ( i = ( n / 2 ) - 1; 0 <= i; i-- )
  {




    key = indx[i];
    ifree = i;

    for ( ; ; )
    {




      m = 2 * ifree + 1;



      if ( n - 1 < m )
      {
        break;
      }



      if ( m + 1 <= n - 1 )
      {




        if ( a[indx[m]] < a[indx[m+1]] )
        {
          m = m + 1;
        }
      }





      if ( a[indx[m]] <= a[key] )
      {
        break;
      }

      indx[ifree] = indx[m];
      ifree = m;
    }



    indx[ifree] = key;
  }

  return;
}


int r8vec_indexed_heap_d_extract ( int &n, double a[], int indx[] )























































{
  int indx_extract;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R8VEC_INDEXED_HEAP_D_EXTRACT - Fatal error!\n";
    cerr << "  The heap is empty.\n";
    exit ( 1 );
  }



  indx_extract = indx[0];

  if ( n == 1 )
  {
    n = 0;
    return indx_extract;
  }



  indx[0] = indx[n-1];



  n = n - 1;
  r8vec_indexed_heap_d ( n, a, indx );

  return indx_extract;
}


void r8vec_indexed_heap_d_insert ( int &n, double a[], int indx[],
  int indx_insert )



















































{
  int i;
  int parent;

  n = n + 1;
  i = n - 1;

  while ( 0 < i )
  {
    parent = ( i - 1 ) / 2;

    if ( a[indx_insert] <= a[indx[parent]] )
    {
      break;
    }

    indx[i] = indx[parent];
    i = parent;
  }

  indx[i] = indx_insert;

  return;
}


int r8vec_indexed_heap_d_max ( int n, double a[], int indx[] )















































{
  int indx_max;

  indx_max = indx[0];

  return indx_max;
}


void r8vec_indicator0 ( int n, double a[] )

























{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = ( double ) ( i );
  }

  return;
}


double *r8vec_indicator0_new ( int n )





























{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = ( double ) ( i );
  }

  return a;
}


void r8vec_indicator1 ( int n, double a[] )

























{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = ( double ) ( i + 1 );
  }

  return;
}


double *r8vec_indicator1_new ( int n )





























{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = ( double ) ( i + 1 );
  }

  return a;
}


void r8vec_insert ( int n, double a[], int pos, double value )




































{
  int i;

  if ( pos < 1 || n + 1 < pos )
  {
    cerr << "\n";
    cerr << "R8VEC_INSERT - Fatal error!\n";
    cerr << "  Illegal insertion position = " << pos << "\n";;
    exit ( 1 );
  }
  else
  {
    for ( i = n + 1; pos + 1 <= i; i-- )
    {
      a[i-1] = a[i-2];
    }

    a[pos-1] = value;
  }

  return;
}


bool r8vec_is_ascending ( int n, double x[] )










































{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < n - 1; i++ )
  {
    if ( x[i+1] < x[i] )
    {
      value = false;
      break;
    }
  }

  return value;
}


bool r8vec_is_ascending_strictly ( int n, double x[] )










































{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < n - 1; i++ )
  {
    if ( x[i+1] <= x[i] )
    {
      value = false;
      break;
    }
  }

  return value;
}


bool r8vec_is_binary ( int n, double x[] )































{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] != 0.0 && x[i] != 1.0 )
    {
      value = false;
      break;
    }
  }
  return value;
}


bool r8vec_is_distinct ( int n, double x[] )































{
  int i;
  int j;
  bool value;

  value = true;

  for ( i = 1; i < n; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      if ( x[i] == x[j] )
      {
        value = false;
        break;
      }
    }
  }
  return value;
}


bool r8vec_is_in_01 ( int n, double x[] )
































{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] < 0.0 || 1.0 < x[i] )
    {
      value = false;
      break;
    }
  }
  return value;
}


bool r8vec_is_in_ab ( int n, double x[], double a, double b )


































{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] < a || b < x[i] )
    {
      value = false;
      break;
    }
  }
  return value;
}


bool r8vec_is_insignificant ( int n, double r[], double s[] )






























{
  int i;
  double t;
  double tol;
  bool value;

  value = true;

  for ( i = 0; i < n; i++ )
  {
    t = r[i] + s[i];
    tol = r8_epsilon ( ) * fabs ( r[i] );

    if ( tol < fabs ( r[i] - t ) )
    {
      value = false;
      break;
    }
  }
  
  return value;
}


bool r8vec_is_integer ( int n, double a[] )































{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] != ( double ) ( int ) a[i] )
    {
      value = false;
      break;
    }
  }
  return value;
}


bool r8vec_is_negative ( int n, double a[] )
































{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < n; i++ )
  {
    if ( 0 <= a[i] )
    {
      value = false;
      break;
    }
  }
  return value;
}


bool r8vec_is_negative_any ( int n, double a[] )
































{
  int i;
  bool value;

  value = false;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < 0.0 )
    {
      value = true;
      break;
    }
  }

  return value;
}


bool r8vec_is_nonnegative ( int n, double x[] )































{
  int i;
  bool value;

  value = true;
  for ( i = 0; i < n; i++ )
  {
    if ( x[i] < 0.0 )
    {
      value = false;
      break;
    }
  }
  return value;
}


bool r8vec_is_nonpositive ( int n, double a[] )
































{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < n; i++ )
  {
    if ( 0.0 < a[i] )
    {
      value = false;
      break;
    }
  }

  return value;
}


bool r8vec_is_nonzero_any ( int n, double a[] )































{
  int i;
  bool value;

  value = false;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] != 0.0 )
    {
      value = true;
      break;
    }
  }

  return value;
}


bool r8vec_is_one ( int n, double x[] )































{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] != 1.0 )
    {
      value = false;
      break;
    }
  }
  return value;
}


bool r8vec_is_positive ( int n, double a[] )
































{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] <= 0.0 )
    {
      value = false;
      break;
    }
  }
  return value;
}


bool r8vec_is_zero ( int n, double x[] )































{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] != 0.0 )
    {
      value = false;
      break;
    }
  }
  return value;
}


double *r8vec_legendre_new ( int n, double a_first, double a_last )































{
  double *a;
  int i;

  a = legendre_zeros ( n );

  for ( i = 0; i < n; i++ )
  {
    a[i] = ( ( 1.0 - a[i] ) * a_first  
           + ( 1.0 + a[i] ) * a_last ) 
           /   2.0;
  }
  return a;
}


void r8vec_linspace ( int n, double a_first, double a_last, double a[] )




































{
  int i;

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - 1 - i ) * a_first 
             + ( double ) (         i ) * a_last ) 
             / ( double ) ( n - 1     );
    }
  }
  return;
}


double *r8vec_linspace_new ( int n, double a_first, double a_last )




































{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - 1 - i ) * a_first 
             + ( double ) (         i ) * a_last ) 
             / ( double ) ( n - 1     );
    }
  }
  return a;
}


double *r8vec_linspace2_new ( int n, double a_first, double a_last )




































{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - i     ) * a_first 
             + ( double ) (     i + 1 ) * a_last ) 
             / ( double ) ( n     + 1 );
    }
  }
  return a;
}


bool r8vec_lt ( int n, double a1[], double a2[] )






































{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a1[i] < a2[i] )
    {
      return true;
    }
    else if ( a2[i] < a1[i] )
    {
      return false;
    }

  }

  return false;
}


void r8vec_mask_print ( int n, double a[], int mask_num, int mask[],
  string title )



































{
  int i;

  cout << "\n";
  cout << "  Masked vector printout:\n";

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < mask_num; i++ )
  {
    cout << "  " << setw(6)  << i
         << ": " << setw(6)  << mask[i]
         << "  " << setw(12) << a[mask[i]-1] << "\n";
  }

  return;
}


double r8vec_max ( int n, double r8vec[] )
































{
  int i;
  double value;

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
  return value;
}


int r8vec_max_abs_index ( int n, double a[] )
































{
  int i;
  int max_index;

  if ( n <= 0 )
  {
    max_index = -1;
  }
  else
  {
    max_index = 0;

    for ( i = 1; i < n; i++ )
    {
      if ( fabs ( a[max_index] ) < fabs ( a[i] ) )
      {
        max_index = i;
      }
    }
  }

  return max_index;
}


int r8vec_max_index ( int n, double a[] )































{
  int i;
  int max_index;

  if ( n <= 0 )
  {
    max_index = -1;
  }
  else
  {
    max_index = 0;

    for ( i = 1; i < n; i++ )
    {
      if ( a[max_index] < a[i] )
      {
        max_index = i;
      }
    }
  }

  return max_index;
}


double r8vec_mean ( int n, double x[] )































{
  int i;
  double mean;

  mean = 0.0;
  for ( i = 0; i < n; i++ )
  {
    mean = mean + x[i];
  }

  mean = mean / ( double ) n;

  return mean;
}


double r8vec_mean_geometric ( int n, double x[] )

































{
  int i;
  double mean;

  mean = 0.0;
  for ( i = 0; i < n; i++ )
  {
    mean = mean + log ( x[i] );
  }

  mean = mean / ( double ) n;
  mean = exp ( mean );

  return mean;
}


double *r8vec_mean_running ( int n, double v[] )




























{
  double *a;
  int i;

  a = new double[n+1];



  a[0] = 0.0;
  for ( i = 1; i < n + 1; i++ )
  {
    a[i] = a[i-1] + v[i-1];
  }



  for ( i = 1; i < n + 1; i++ )
  {
    a[i] = a[i] / ( double ) ( i );
  }

  return a;
}


double r8vec_median ( int n, double a[] )



































{
  int k;
  double median;

  k = ( n + 1 ) / 2;

  median = r8vec_frac ( n, a, k );

  return median;
}


void r8vec_mesh_2d ( int nx, int ny, double xvec[], double yvec[], 
  double xmat[], double ymat[] )










































{
  int i;
  int j;

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      xmat[i+j*nx] = xvec[i];
    }
  }

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      ymat[i+j*nx] = yvec[j];
    }
  }

 return;
}


double *r8vec_midspace_new ( int n, double a, double b )







































{
  double *x;
  int i;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( 2 * n - 2 * i - 1 ) * a 
           + ( double ) (         2 * i + 1 ) * b ) 
           / ( double ) ( 2 * n );
  }

  return x;
}


double r8vec_min ( int n, double r8vec[] )































{
  int i;
  double value;

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
  }
  return value;
}


int r8vec_min_index ( int n, double a[] )































{
  int i;
  int min_index;

  if ( n <= 0 )
  {
    min_index = -1;
  }
  else
  {
    min_index = 0;

    for ( i = 1; i < n; i++ )
    {
      if ( a[i] < a[min_index] )
      {
        min_index = i;
      }
    }
  }

  return min_index;
}


double r8vec_min_pos ( int n, double a[] )































{
  int i;
  const double r8_huge = 1.79769313486231571E+308;
  double value;

  value = r8_huge;

  for ( i = 0; i < n; i++ )
  {
    if ( 0.0 < a[i] )
    {
      if ( a[i] < value )
      {
        value = a[i];
      }
    }
  }
  return value;
}


bool r8vec_mirror_next ( int n, double a[] )























































































{
  bool done;
  int i;
  int positive;



  positive = -1;
  for ( i = 0; i < n; i++ )
  {
    if ( 0.0 < a[i] )
    {
      positive = i;
      break;
    }
  }



  if ( positive == -1 )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = - a[i];
    }
    done = true;
    return done;
  }



  for ( i = 0; i <= positive; i++ )
  {
    a[i] = - a[i];
  }
  done = false;

  return done;
}


void r8vec_mirror_ab_next ( int m, double a[], double b[], double x[], 
  bool &done )









































































































{
  int i;



  if ( done )
  {



    for ( i = 0; i < m; i++ )
    {
      if ( x[i] < a[i] )
      {
        cerr << "\n";
        cerr << "R8VEC_MIRROR_AB_NEXT - Fatal error!\n";
        cerr << "  Not every A(I) <= X(I).\n";
        exit ( 1 );
      }
      if ( b[i] < x[i] )
      {
        cerr << "\n";
        cerr << "R8VEC_MIRROR_AB_NEXT - Fatal error!\n";
        cerr << "  Not every X(I) <= B(I).\n";
        exit ( 1 );
      }
    }



    for ( i = 0; i < m; i++ )
    {
      x[i] = 2.0 * a[i] - x[i];
    }



    done = true;
    for ( i = 0; i < m; i++ )
    {
      if ( a[i] != b[i] )
      {
        done = false;
        break;
      }
    }
  }



  else
  {








    i = m - 1;

    while ( 0 <= i )
    {
      if ( x[i] < a[i] )
      {
        x[i] = 2.0 * a[i] - x[i];
        return;
      }
      else if ( x[i] < b[i] )
      {
        x[i] = 2.0 * b[i] - x[i];
        return;
      }
      else
      {
        x[i] = x[i] - 2.0 * ( b[i] - a[i] );
      }
      i = i - 1;
    }
    done = true;
  }

  return;
}


void r8vec_mm_to_01 ( int n, double a[] )




































{
  double amax;
  double amin;
  int i;

  amax = r8vec_max ( n, a );
  amin = r8vec_min ( n, a );

  if ( amin == amax )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = 0.5;
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( a[i] - amin ) / ( amax - amin );
    }
  }

  return;
}


double *r8vec_mm_to_cd ( int n, double a[], double bmin, double bmax )





































{
  double amax;
  double amin;
  double *b;
  int i;

  b = new double[n];

  if ( bmax == bmin )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i] = bmin;
    }
    return b;
  }

  amax = r8vec_max ( n, a );
  amin = r8vec_min ( n, a );

  if ( amin == amax )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i] = 0.5 * ( bmax + bmin );
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      b[i] = ( ( amax - a[i]        ) * bmin
             + (        a[i] - amin ) * bmax )
             / ( amax        - amin );
    }
  }

  return b;
}


void r8vec_nint ( int n, double a[] )





























{
  int i;
  int s;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < 0.0 )
    {
      s = -1;
    }
    else
    {
      s = 1;
    }
    a[i] = ( double ) ( s * ( int ) ( fabs ( a[i] ) + 0.5 ) );
  }

  return;
}


double *r8vec_nint_new ( int n, double a[] )































{
  double *b;
  int i;
  int s;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < 0.0 )
    {
      s = -1;
    }
    else
    {
      s = 1;
    }
    b[i] = ( double ) ( s * ( int ) ( fabs ( a[i] ) + 0.5 ) );
  }

  return b;
}


double r8vec_norm ( int n, double a[] )



































{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  v = sqrt ( v );

  return v;
}


double r8vec_norm_affine ( int n, double v0[], double v1[] )




































{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( v1[i] - v0[i] ) * ( v1[i] - v0[i] );
  }
  value = sqrt ( value );

  return value;
}


double r8vec_norm_l0 ( int n, double a[] )



































{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    if ( a[i] != 0.0 )
    {
      value = value + 1.0;
    }
  }
  return value;
}


double r8vec_norm_l1 ( int n, double a[] )



































{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + fabs ( a[i] );
  }

  return v;
}


double r8vec_norm_l2 ( int n, double a[] )



































{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  v = sqrt ( v );

  return v;
}


double r8vec_norm_li ( int n, double a[] )



































{
  int i;
  double v1;
  double v2;

  v1 = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v2 = fabs ( a[i] );

    if ( v1 < v2 )
    {
      v1 = v2;
    }
  }

  return v1;
}


double r8vec_norm_lp ( int n, double a[], double p )





































{
  int i;
  double v;

  v = 0.0;

  if ( p == 1.0 )
  {
    for ( i = 0; i < n; i++ )
    {
      v = v + fabs ( a[i] );
    }
  }
  else if ( p == 2.0 )
  {
    for ( i = 0; i < n; i++ )
    {
      v = v + a[i] * a[i];
    }
    v = sqrt ( v );
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      v = v + pow ( fabs ( a[i] ), p );
    }
    v = pow (  ( double ) v, 1.0 / p );
  }

  return v;
}


void r8vec_normal_01 ( int n, int &seed, double x[] )











































{
  int i;
  int m;
  double *r;
  const double r8_pi = 3.141592653589793;
  int x_hi;
  int x_lo;



  x_lo = 1;
  x_hi = n;



  if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * r8_pi * r[1] );

    delete [] r;
  }



  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    delete [] r;
  }





  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );

    delete [] r;
  }

  return;
}


double *r8vec_normal_01_new ( int n, int &seed )











































{
  int i;
  int m;
  double *r;
  const double r8_pi = 3.141592653589793;
  double *x;
  int x_hi;
  int x_lo;

  x = new double[n];



  x_lo = 1;
  x_hi = n;



  if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( - 2.0 * log ( r[0] ) ) * cos ( 2.0 * r8_pi * r[1] );

    delete [] r;
  }



  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2 * m, seed );

    for ( i = 0; i <= 2 * m - 2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    delete [] r;
  }





  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2 * m, seed );

    for ( i = 0; i <= 2 * m - 4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    i = 2 * m - 2;

    x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );

    delete [] r;
  }

  return x;
}


double *r8vec_normal_ab_new ( int n, double b, double c, int &seed )













































{
  int i;
  int m;
  double *r;
  const double r8_pi = 3.141592653589793;
  double *x;
  int x_hi;
  int x_lo;

  x = new double[n];



  x_lo = 1;
  x_hi = n;



  if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( - 2.0 * log ( r[0] ) ) * cos ( 2.0 * r8_pi * r[1] );

    delete [] r;
  }



  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2 * m, seed );

    for ( i = 0; i <= 2 * m - 2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    delete [] r;
  }





  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2 * m - 4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );

    delete [] r;
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = b + c * x[i];
  }

  return x;
}


void r8vec_normalize ( int n, double a[] )






























{
  int i;
  double norm;

  norm = 0.0;
  for ( i = 0; i < n; i++ )
  {
    norm = norm + a[i] * a[i];
  }
  norm = sqrt ( norm );

  if ( norm == 0.0 )
  {
    cerr << "\n";
    cerr << "R8VEC_NORMALIZE - Fatal error!\n";
    cerr << "  The vector norm is 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] / norm;
  }

  return;
}


void r8vec_normalize_l1 ( int n, double a[] )































{
  double a_sum;
  int i;

  a_sum = 0.0;
  for ( i = 0; i < n; i++ )
  {
    a_sum = a_sum + a[i];
  }

  if ( a_sum == 0.0 )
  {
    cerr << "\n";
    cerr << "R8VEC_NORMALIZE_L1 - Fatal error!\n";
    cerr << "  The vector entries sum to 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] / a_sum;
  }

  return;
}


double r8vec_normsq ( int n, double a[] )



































{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  return v;
}


double r8vec_normsq_affine ( int n, double v0[], double v1[] )




































{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( v1[i] - v0[i] ) * ( v1[i] - v0[i] );
  }
  return value;
}


double *r8vec_ones_new ( int n )





























{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 1.0;
  }
  return a;
}


int r8vec_order_type ( int n, double x[] )





































{
  int i;
  int order;



  i = 0;

  for ( ; ; )
  {
    i = i + 1;
    if ( n-1 < i )
    {
      order = 0;
      return order;
    }

    if ( x[0] < x[i] )
    {
      if ( i == 1 )
      {
        order = 2;
        break;
      }
      else
      {
        order = 1;
        break;
      }
    }
    else if ( x[i] < x[0] )
    {
      if ( i == 1 )
      {
        order = 4;
        break;
      }
      else
      {
        order = 3;
        break;
      }
    }
  }



  for ( ; ; )
  {
    i = i + 1;
    if ( n - 1 < i )
    {
      break;
    }

    if ( order == 1 )
    {
      if ( x[i] < x[i-1] )
      {
        order = -1;
        break;
      }
    }
    else if ( order == 2 )
    {
      if ( x[i] < x[i-1] )
      {
        order = -1;
        break;
      }
      else if ( x[i] == x[i-1] )
      {
        order = 1;
      }
    }
    else if ( order == 3 )
    {
      if ( x[i-1] < x[i] )
      {
        order = -1;
        break;
      }
    }
    else if ( order == 4 )
    {
      if ( x[i-1] < x[i] )
      {
        order = -1;
        break;
      }
      else if ( x[i] == x[i-1] )
      {
        order = 3;
      }
    }
  }
  return order;
}


void r8vec_part_quick_a ( int n, double a[], int &l, int &r )























































{
  int i;
  double key;
  int m;
  double temp;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R8VEC_PART_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }
  else if ( n == 1 )
  {
    l = 0;
    r = 2;
    return;
  }

  key = a[0];
  m = 1;



  l = 1;
  r = n + 1;

  for ( i = 2; i <= n; i++ )
  {

    if ( key < a[l] )
    {
      r = r - 1;
      temp = a[r-1];
      a[r-1] = a[l];
      a[l] = temp;
    }
    else if ( a[l] == key )
    {
      m = m + 1;
      temp = a[m-1];
      a[m-1] = a[l];
      a[l] = temp;
      l = l + 1;
    }
    else if ( a[l] < key )
    {
      l = l + 1;
    }

  }



  for ( i = 1; i <= l - m; i++ )
  {
    a[i-1] = a[i+m-1];
  }

  l = l - m;

  for ( i = l + 1; i <= l + m; i++ )
  {
    a[i-1] = key;
  }

  return;
}


void r8vec_permute ( int n, int p[], double a[] )


















































{
  double a_temp;
  int i;
  int iget;
  int iput;
  int istart;

  if ( !perm0_check ( n, p ) )
  {
    cerr << "\n";
    cerr << "R8VEC_PERMUTE - Fatal error!\n";
    cerr << "  PERM0_CHECK rejects permutation.\n";
    exit ( 1 );
  }





  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1;
  }



  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      a_temp = a[istart-1];
      iget = istart;



      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cerr << "\n";
          cerr << "R8VEC_PERMUTE - Fatal error!\n";
          cerr << "  A permutation index is out of range.\n";
          cerr << "  P(" << iput << ") = " << iget << "\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[iput-1] = a_temp;
          break;
        }
        a[iput-1] = a[iget-1];
      }
    }
  }



  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }



  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1;
  }
  return;
}


void r8vec_permute_cyclic ( int n, int k, double a[] )




































{
  double *b;
  int i;
  int ipk;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    ipk = i4_wrap ( i + k, 0, n - 1 );
    b[i] = a[ipk];
  }

  for ( i = 0; i < n; i++ )
  {
    a[i] = b[i];
  }

  delete [] b;

  return;
}


void r8vec_permute_uniform ( int n, double a[], int &seed )































{
  int *p;

  p = perm0_uniform_new ( n, seed );

  r8vec_permute ( n, p, a );

  delete [] p;

  return;
}


void r8vec_polarize ( int n, double a[], double p[], double a_normal[],
  double a_parallel[] )











































{
  double a_dot_p;
  int i;
  double p_norm;

  p_norm = r8vec_norm ( n, p );

  if ( p_norm == 0.0 )
  {
    for ( i = 0; i < n; i++ )
    {
      a_normal[i] = a[i];
    }
    for ( i = 0; i < n; i++ )
    {
      a_parallel[i] = 0.0;
    }
    return;
  }
  a_dot_p = r8vec_dot_product ( n, a, p ) / p_norm;

  for ( i = 0; i < n; i++ )
  {
    a_parallel[i] = a_dot_p * p[i] / p_norm;
  }

  for ( i = 0; i < n; i++ )
  {
    a_normal[i] = a[i] - a_parallel[i];
  }

  return;
}


void r8vec_print ( int n, double a[], string title )































{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8)  << i
         << ": " << setw(14) << a[i]  << "\n";
  }

  return;
}


void r8vec_print_16 ( int n, double a[], string title )































{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8)  << i
         << ": " << setprecision(16) << setw(24) << a[i]  << "\n";
  }

  return;
}


void r8vec_print_part ( int n, double a[], int i_lo, int i_hi, string title )


































{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = i4_max ( 1, i_lo ); i <= i4_min ( n, i_hi ); i++ )
  {
    cout << "  " << setw(8)  << i
         << ": " << setw(14) << a[i-1]  << "\n";
  }

  return;
}


void r8vec_print_some ( int n, double a[], int max_print, string title )









































{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(14) << a[i] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[i] << "\n";
    }
    cout << "  ........  ..............\n";
    i = n - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[i] << "\n";
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[i] << "\n";
    }
    i = max_print - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[i]
         << "  " << "...more entries...\n";
  }

  return;
}


double r8vec_product ( int n, double a[] )































{
  int i;
  double product;

  product = 1.0;
  for ( i = 0; i < n; i++ )
  {
    product = product * a[i];
  }

  return product;
}


void r8vec_range ( int n, double x[], double xmin, double xmax, double y[],
  double *ymin, double *ymax )









































{
  int i;
  const double r8_huge = 1.79769313486231571E+308;

  *ymin =   r8_huge;
  *ymax = - r8_huge;

  for ( i = 0; i < n; i++ )
  {
    if ( xmin <= x[i] && x[i] <= xmax )
    {
      *ymin = r8_min ( *ymin, y[i] );
      *ymax = r8_max ( *ymax, y[i] );
    }
  }

  return;
}


void r8vec_range_2 ( int n, double a[], double *amin, double *amax )











































{
  *amax = r8_max ( *amax, r8vec_max ( n, a ) );
  *amin = r8_min ( *amin, r8vec_min ( n, a ) );

  return;
}


void r8vec_reverse ( int n, double a[] )







































{
  int i;
  int i_hi;
  double temp;

  i_hi = n / 2;

  for ( i = 1; i <= i_hi; i++ )
  {
    temp   = a[i-1];
    a[i-1] = a[n-i];
    a[n-i] = temp;
  }

  return;
}


double r8vec_rms ( int n, double a[] )



































{
  int i;
  double v;

  v = 0.0;

  if ( 0 < n )
  {
    for ( i = 0; i < n; i++ )
    {
      v = v + a[i] * a[i];
    }
    v = sqrt ( v / ( double ) ( n ) );
  }
  return v;
}


void r8vec_rotate ( int n, double a[], int m )



















































{
  int iget;
  int iput;
  int istart;
  int mcopy;
  int nset;
  double temp;



  mcopy = i4_modp ( m, n );

  if ( mcopy == 0 )
  {
    return;
  }

  istart = 0;
  nset = 0;

  for ( ; ; )
  {
    istart = istart + 1;

    if ( n < istart )
    {
      break;
    }

    temp = a[istart-1];
    iget = istart;



    for ( ; ; )
    {
      iput = iget;

      iget = iget - mcopy;
      if ( iget < 1 )
      {
        iget = iget + n;
      }

      if ( iget == istart )
      {
        break;
      }

      a[iput-1] = a[iget-1];
      nset = nset + 1;
    }

    a[iput-1] = temp;
    nset = nset + 1;

    if ( n <= nset )
    {
      break;
    }
  }

  return;
}


double r8vec_scalar_triple_product ( double v1[3], double v2[3], double v3[3] )

































{
  double value;

  value =
      v1[0] * ( v2[1] * v3[2] - v2[2] * v3[1] )
    + v1[1] * ( v2[2] * v3[0] - v2[0] * v3[2] )
    + v1[2] * ( v2[0] * v3[1] - v2[1] * v3[0] );

  return value;
}


void r8vec_scale ( double s, int n, double a[] )
































{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = s * a[i];
  }
  return;
}


int r8vec_search_binary_a ( int n, double a[], double aval )













































{
  int high;
  int indx;
  int low;
  int mid;

  indx = -1;

  low = 1;
  high = n;

  while ( low <= high )
  {
    mid = ( low + high ) / 2;

    if ( a[mid-1] == aval )
    {
      indx = mid;
      break;
    }
    else if ( a[mid-1] < aval )
    {
      low = mid + 1;
    }
    else if ( aval < a[mid-1] )
    {
      high = mid - 1;
    }
  }

  return indx;
}


void r8vec_shift ( int shift, int n, double x[] )
































{
  int i;
  int ihi;
  int ilo;
  double *y;

  y = new double[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = x[i];
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }

  ilo = i4_max ( 0, shift );
  ihi = i4_min ( n, n + shift );

  for ( i = ilo; i < ihi; i++ )
  {
    x[i] = y[i-shift];
  }

  delete [] y;

  return;
}


void r8vec_shift_circular ( int shift, int n, double x[] )
































{
  int i;
  int j;
  double *y;

  y = new double[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = x[i];
  }

  for ( i = 0; i < n; i++ )
  {
    j = i4_wrap ( i - shift, 0, n - 1 );
    x[i] = y[j];
  }
  delete [] y;
  return;
}


double *r8vec_sign3_running ( int n, double v[] )































{
  int i;
  double *s;

  s = new double[n+1];



  s[0] = 0.0;
  for ( i = 1; i < n + 1; i++ )
  {
    s[i] = s[i-1] + v[i-1];
  }

  for ( i = 0; i < n + 1; i++ )
  {
    if ( s[i] < 0.0 )
    {
      s[i] = -1.0;
    }
    else if ( s[i] == 0.0 )
    {
      s[i] = 0.0;
    }
    else if ( 0.0 < s[i] )
    {
      s[i] = +1.0;
    }
  }

  return s;
}


void r8vec_sort_bubble_a ( int n, double a[] )































{
  int i;
  int j;
  double temp;

  for ( i = 0; i < n-1; i++ )
  {
    for ( j = i+1; j < n; j++ )
    {
      if ( a[j] < a[i] )
      {
        temp = a[i];
        a[i] = a[j];
        a[j] = temp;
      }
    }
  }
  return;
}


void r8vec_sort_bubble_d ( int n, double a[] )































{
  int i;
  int j;
  double temp;

  for ( i = 0; i < n-1; i++ )
  {
    for ( j = i+1; j < n; j++ )
    {
      if ( a[i] < a[j] )
      {
        temp = a[i];
        a[i] = a[j];
        a[j] = temp;
      }
    }
  }
  return;
}


void r8vec_sort_heap_a ( int n, double a[] )






































{
  int n1;
  double temp;

  if ( n <= 1 )
  {
    return;
  }



  r8vec_heap_d ( n, a );






  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;



  for ( n1 = n-1; 2 <= n1; n1-- )
  {



    r8vec_heap_d ( n1, a );



    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;
  }

  return;
}


void r8vec_sort_heap_d ( int n, double a[] )






































{
  int n1;
  double temp;

  if ( n <= 1 )
  {
    return;
  }



  r8vec_heap_a ( n, a );






  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;



  for ( n1 = n-1; 2 <= n1; n1-- )
  {



    r8vec_heap_a ( n1, a );



    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;
  }

  return;
}


void r8vec_sort_heap_index_a ( int n, double a[], int indx[] )


















































{
  double aval;
  int i;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return;
  }

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    return;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j-1]] < a[indx[j]] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[indx[j-1]] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return;
}


int *r8vec_sort_heap_index_a_new ( int n, double a[] )


















































{
  double aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = new int[n];

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j-1]] < a[indx[j]] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[indx[j-1]] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return indx;
}


void r8vec_sort_heap_index_d ( int n, double a[], int indx[] )


















































{
  double aval;
  int i;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return;
  }

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    return;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j]] < a[indx[j-1]] )
        {
          j = j + 1;
        }
      }

      if ( a[indx[j-1]] < aval )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }

    indx[i-1] = indxt;
  }
  return;
}


int *r8vec_sort_heap_index_d_new ( int n, double a[] )


















































{
  double aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = new int[n];

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j]] < a[indx[j-1]] )
        {
          j = j + 1;
        }
      }

      if ( a[indx[j-1]] < aval )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }

    indx[i-1] = indxt;
  }
  return indx;
}


int *r8vec_sort_heap_mask_a ( int n, double a[], int mask_num, int mask[] )






















































{
  double aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  if ( mask_num < 1 )
  {
    return NULL;
  }

  if ( mask_num == 1 )
  {
    indx = new int[1];
    indx[0] = 0;
    return indx;
  }

  indx = i4vec_indicator1_new ( mask_num );

  l = mask_num / 2 + 1;
  ir = mask_num;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[mask[indxt-1]-1];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[mask[indxt-1]-1];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[mask[indx[j-1]-1]-1] < a[mask[indx[j]-1]-1] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[mask[indx[j-1]-1]-1] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  for ( i = 0; i < mask_num; i++ )
  {
    indx[i] = indx[i] - 1;
  }

  return indx;
}


void r8vec_sort_insert_a ( int n, double a[] )








































{
  int i;
  int j;
  double x;

  for ( i = 1; i < n; i++ )
  {
    x = a[i];

    j = i;

    while ( 1 <= j && x < a[j-1] )
    {
      a[j] = a[j-1];
      j = j - 1;
    }
    a[j] = x;
  }

  return;
}


int *r8vec_sort_insert_index_a ( int n, double a[] )







































{
  int i;
  int *indx;
  int j;
  double x;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = i4vec_indicator0_new ( n );

  for ( i = 1; i < n; i++ )
  {
    x = a[i];
    j = i - 1;

    while ( 0 <= j )
    {
      if ( a[indx[j]] <= x )
      {
        break;
      }

      indx[j+1] = indx[j];
      j = j - 1;
    }
    indx[j+1] = i;
  }

  return indx;
}


void r8vec_sort_quick_a ( int n, double a[] )










































{
# define LEVEL_MAX 30

  int base;
  int l_segment;
  int level;
  int n_segment;
  int rsave[LEVEL_MAX];
  int r_segment;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R8VEC_SORT_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }
  else if ( n == 1 )
  {
    return;
  }

  level = 1;
  rsave[0] = n + 1;
  base = 1;
  n_segment = n;

  while ( 0 < n_segment )
  {



    r8vec_part_quick_a ( n_segment, a+base-1, l_segment, r_segment );



    if ( 1 < l_segment )
    {

      if ( LEVEL_MAX < level )
      {
        cerr << "\n";
        cerr << "R8VEC_SORT_QUICK_A - Fatal error!\n";
        cerr << "  Exceeding recursion maximum of " << LEVEL_MAX << "\n";
        exit ( 1 );
      }

      level = level + 1;
      n_segment = l_segment;
      rsave[level-1] = r_segment + base - 1;
    }




    else if ( r_segment < n_segment )
    {
      n_segment = n_segment + 1 - r_segment;
      base = base + r_segment - 1;
    }



    else
    {
      for ( ; ; )
      {
        if ( 1 < level )
        {
          base = rsave[level-1];
          n_segment = rsave[level-2] - rsave[level-1];
          level = level - 1;
          if ( 0 < n_segment )
          {
            break;
          }
        }
        else
        {
          n_segment = 0;
          break;
        }
      }
    }
  }

  return;
# undef LEVEL_MAX
}


void r8vec_sort_shell_a ( int n, double a[] )































{
  double asave;
  int i;
  int ifree;
  int inc;
  int ipow;
  int j;
  int k;
  int maxpow;
  int test;

  if ( n <= 1 )
  {
    return;
  }




  maxpow = 1;
  test = 3;

  while ( test < 2 * n + 1 )
  {
    maxpow = maxpow + 1;
    test = test * 3;
  }

  if ( 1 < maxpow )
  {
    maxpow = maxpow - 1;
    test = test / 3;
  }



  for ( ipow = maxpow; 1 <= ipow; ipow-- )
  {
    inc = ( test - 1 ) / 2;
    test = test / 3;



    for ( k = 1; k <= inc; k++ )
    {




      for ( i = inc+k; i <= n; i = i + inc )
      {
        asave = a[i-1];
        ifree = i;
        j = i - inc;

        for ( ; ; )
        {
          if ( j < 1 )
          {
            break;
          }

          if ( a[j-1] <= asave )
          {
            break;
          }

          ifree = j;
          a[j+inc-1] = a[j-1];
          j = j - inc;
        }
        a[ifree-1] = asave;
      }
    }
  }

  return;
}


double *r8vec_sorted_merge_a ( int na, double a[], int nb, double b[], int &nc )










































{
  double *c;
  double *d;
  int j;
  int ja;
  int jb;
  int na2;
  int nb2;
  int nd;
  int order;

  na2 = na;
  nb2 = nb;

  ja = 0;
  jb = 0;
  nc = 0;
  nd = 0;
  d = new double[na+nb];

  order = r8vec_order_type ( na2, a );

  if ( order < 0 || 2 < order )
  {
    cerr << "\n";
    cerr << "R8VEC_SORTED_MERGE_A - Fatal error!\n";
    cerr << "  The input array A is not ascending sorted.\n";
    return NULL;
  }

  order = r8vec_order_type ( nb2, b );

  if ( order < 0 || 2 < order )
  {
    cerr << "\n";
    cerr << "R8VEC_SORTED_MERGE_A - Fatal error!\n";
    cerr << "  The input array B is not ascending sorted.\n";
    return NULL;
  }

  for ( ; ; )
  {



    if ( na2 <= ja )
    {
      for ( j = 1; j <= nb2 - jb; j++ )
      {
        jb = jb + 1;
        if ( nd == 0 )
        {
          nd = nd + 1;
          d[nd-1] = b[jb-1];
        }
        else if ( d[nd-1] < b[jb-1] )
        {
          nd = nd + 1;
          d[nd-1] = b[jb-1];
        }
      }
      break;
    }



    else if ( nb2 <= jb )
    {
      for ( j = 1; j <= na2 - ja; j++ )
      {
        ja = ja + 1;
        if ( nd == 0 )
        {
          nd = nd + 1;
          d[nd-1] = a[ja-1];
        }
        else if ( d[nd-1] < a[ja-1] )
        {
          nd = nd + 1;
          d[nd-1] = a[ja-1];
        }
      }
      break;
    }



    else if ( a[ja] <= b[jb] )
    {
      ja = ja + 1;
      if ( nd == 0 )
      {
        nd = nd + 1;
        d[nd-1] = a[ja-1];
      }
      else if ( d[nd-1] < a[ja-1] )
      {
        nd = nd + 1;
        d[nd-1] = a[ja-1];
      }
    }



    else
    {
      jb = jb + 1;
      if ( nd == 0 )
      {
        nd = nd + 1;
        d[nd-1] = b[jb-1];
      }
      else if ( d[nd-1] < b[jb-1] )
      {
        nd = nd + 1;
        d[nd-1] = b[jb-1];
      }
    }
  }

  nc = nd;

  c = r8vec_copy_new ( nd, d );

  delete [] d;

  return c;
}


int r8vec_sorted_nearest ( int n, double a[], double value )


































{
  int hi;
  int lo;
  int mid;

  if ( n < 1 )
  {
    return (-1);
  }

  if ( n == 1 )
  {
    return 1;
  }

  if ( a[0] < a[n-1] )
  {
    if ( value < a[0] )
    {
      return 1;
    }
    else if ( a[n-1] < value )
    {
      return n;
    }



    lo = 1;
    hi = n;

    while ( lo < hi - 1 )
    {
      mid = ( lo + hi ) / 2;

      if ( value == a[mid-1] )
      {
        return mid;
      }
      else if ( value < a[mid-1] )
      {
        hi = mid;
      }
      else
      {
        lo = mid;
      }
    }



    if ( fabs ( value - a[lo-1] ) < fabs ( value - a[hi-1] ) )
    {
      return lo;
    }
    else
    {
      return hi;
    }
  }



  else
  {
    if ( value < a[n-1] )
    {
      return n;
    }
    else if ( a[0] < value )
    {
      return 1;
    }



    lo = n;
    hi = 1;

    while ( lo < hi - 1 )
    {
      mid = ( lo + hi ) / 2;

      if ( value == a[mid-1] )
      {
        return mid;
      }
      else if ( value < a[mid-1] )
      {
        hi = mid;
      }
      else
      {
        lo = mid;
      }
    }



    if ( fabs ( value - a[lo-1] ) < fabs ( value - a[hi-1] ) )
    {
      return lo;
    }
    else
    {
      return hi;
    }
  }
}


void r8vec_sorted_range ( int n, double r[], double r_lo, double r_hi,
  int &i_lo, int &i_hi )































{
  int i1;
  int i2;
  int j1;
  int j2;



  if ( r[n-1] < r_lo )
  {
    i_lo = - 1;
    i_hi = - 2;
    return;
  }

  if ( r_hi < r[0] )
  {
    i_lo = - 1;
    i_hi = - 2;
    return;
  }



  if ( n == 1 )
  {
    if ( r_lo <= r[0] && r[0] <= r_hi )
    {
      i_lo = 1;
      i_hi = 1;
    }
    else
    {
      i_lo = - 1;
      i_hi = - 2;
    }
    return;
  }



  if ( r_lo <= r[0] )
  {
    i_lo = 0;
  }
  else
  {





    j1 = 0;
    j2 = n - 1;
    i1 = ( j1 + j2 - 1 ) / 2;
    i2 = i1 + 1;

    for ( ; ; )
    {
      if ( r_lo < r[i1] )
      {
        j2 = i1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else if ( r[i2] < r_lo )
      {
        j1 = i2;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else
      {
        i_lo = i1;
        break;
      }
    }
  }



  if ( r[n-1] <= r_hi )
  {
    i_hi = n - 1;
  }
  else
  {
    j1 = i_lo;
    j2 = n - 1;
    i1 = ( j1 + j2 - 1 ) / 2;
    i2 = i1 + 1;

    for ( ; ; )
    {
      if ( r_hi < r[i1] )
      {
        j2 = i1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else if ( r[i2] < r_hi )
      {
        j1 = i2;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else
      {
        i_hi = i2;
        break;
      }
    }
  }







  if ( r[i_lo] < r_lo )
  {
    i_lo = i_lo + 1;
    if ( n - 1 < i_lo )
    {
      i_hi = i_lo - 1;
    }
  }

  if ( r_hi < r[i_hi] )
  {
    i_hi = i_hi - 1;
    if ( i_hi < 0 )
    {
      i_lo = i_hi + 1;
    }
  }

  return;
}


void r8vec_sorted_split ( int n, double a[], double split, int &i_lt,
  int &i_gt )


















































{
  int hi;
  int i;
  int lo;
  int mid;

  if ( n < 1 )
  {
    i_lt = -1;
    i_gt = -1;
    return;
  }

  if ( split < a[0] )
  {
    i_lt = 0;
    i_gt = 1;
    return;
  }

  if ( a[n-1] < split )
  {
    i_lt = n;
    i_gt = n + 1;
    return;
  }

  lo = 1;
  hi = n;

  for ( ; ; )
  {
    if ( lo + 1 == hi )
    {
      i_lt = lo;
      break;
    }

    mid = ( lo + hi ) / 2;

    if ( split <= a[mid-1] )
    {
      hi = mid;
    }
    else
    {
      lo = mid;
    }
  }

  for ( i = i_lt + 1; i <= n; i++ )
  {
    if ( split < a[i-1] )
    {
      i_gt = i;
      return;
    }
  }

  i_gt = n + 1;

  return;
}


void r8vec_sorted_undex ( int x_num, double x_val[], int x_unique_num,
  double tol, int undx[], int xdnu[] )




























































































{
  int i;
  int j;



  i = 0;

  j = 0;
  undx[j] = i;

  xdnu[i] = j;

  for ( i = 1; i < x_num; i++ )
  {
    if ( tol < fabs ( x_val[i] - x_val[undx[j]] ) )
    {
      j = j + 1;
      undx[j] = i;
    }
    xdnu[i] = j;
  }

  return;
}


double *r8vec_sorted_unique ( int n, double a[], double tol, int &unique_num )






































{
  double *a_unique;
  int i;
  int iuniq;

  unique_num = 0;

  if ( n <= 0 )
  {
    return NULL;
  }



  iuniq = 0;
  unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( tol < fabs ( a[i] - a[iuniq] ) )
    {
       iuniq = i;
      unique_num = unique_num + 1;
    }
  }



  a_unique = new double[unique_num];



  unique_num = 0;

  a_unique[unique_num] = a[0];
  unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( tol < fabs ( a[i] - a_unique[unique_num-1] ) )
    {
      a_unique[unique_num] = a[i];
      unique_num = unique_num + 1;
    }
  }

  return a_unique;
}


int r8vec_sorted_unique_count ( int n, double a[], double tol )



































{
  int i;
  int unique_num;

  unique_num = 0;

  if ( n < 1 )
  {
    return unique_num;
  }

  unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( tol < fabs ( a[i-1] - a[i] ) )
    {
      unique_num = unique_num + 1;
    }
  }

  return unique_num;
}


void r8vec_sorted_unique_hist ( int n, double a[], double tol, int maxuniq,
  int &unique_num, double auniq[], int acount[] )











































{
  int i;
  int index;



  index = -1;

  for ( i = 0; i < n; i++ )
  {

    if ( i == 0 )
    {
      index = 0;
      auniq[index] = a[0];
      acount[index] = 1;
    }
    else if ( fabs ( a[i] - auniq[index] ) <= tol )
    {
      acount[index] = acount[index] + 1;
    }
    else if ( index + 1 < maxuniq )
    {
      index = index + 1;
      auniq[index] = a[i];
      acount[index] = 1;
    }
  }

  unique_num = index + 1;

  return;
}


int r8vec_split ( int n, double a[], double split )












































{
  int i;
  int i1;
  int i2;
  int i3;
  int isplit;
  int j1;
  int j2;
  int j3;
  double temp;






  i1 = 1;
  j1 = 0;

  i2 = 1;
  j2 = n;

  i3 = n+1;
  j3 = n;




  for ( i = 1; i <= n; i++ )
  {
    if ( a[i2-1] <= split )
    {
      i2 = i2 + 1;
      j1 = j1 + 1;
    }
    else
    {
      temp = a[i2-1];
      a[i2-1] = a[i3-2];
      a[i3-2] = temp;
      i3 = i3 - 1;
      j2 = j2 - 1;
    }
  }

  isplit = j1;

  return isplit;
}


double r8vec_std ( int n, double a[] )






































{
  int i;
  double mean;
  double std;

  if ( n < 2 )
  {
    std = 0.0;
  }
  else
  {
    mean = 0.0;
    for ( i = 0; i < n; i++ )
    {
      mean = mean + a[i];
    }
    mean = mean / ( ( double ) n );

    std = 0.0;
    for ( i = 0; i < n; i++ )
    {
      std = std + ( a[i] - mean ) * ( a[i] - mean );
    }
    std = sqrt ( std / ( ( double ) ( n ) ) );
  }

  return std;
}


double r8vec_std_sample ( int n, double a[] )






































{
  int i;
  double mean;
  double std;

  if ( n < 2 )
  {
    std = 0.0;
  }
  else
  {
    mean = 0.0;
    for ( i = 0; i < n; i++ )
    {
      mean = mean + a[i];
    }
    mean = mean / ( ( double ) n );

    std = 0.0;
    for ( i = 0; i < n; i++ )
    {
      std = std + ( a[i] - mean ) * ( a[i] - mean );
    }
    std = sqrt ( std / ( ( double ) ( n - 1 ) ) );
  }

  return std;
}


void r8vec_std_update ( int nm1, double mean_nm1, double std_nm1, double xn, 
  int &n, double &mean_n, double &std_n )







































{
  if ( nm1 <= 0 )
  {
    n = 1;
    mean_n = xn;
    std_n = 0.0;
  }
  else
  {
    n = nm1 + 1;
    mean_n = mean_nm1 + ( xn - mean_nm1 ) / ( double ) ( n );
    std_n = sqrt ( ( std_nm1 * std_nm1 * ( double ) ( nm1 ) 
      + ( xn - mean_nm1 ) * ( xn - mean_n ) ) / ( double ) ( n ) );
  }

  return;
}


void r8vec_step ( double x0, int n, double x[], double fx[] )


































{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] < x0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = 1.0;
    }
  }
  return;
}


void r8vec_stutter ( int n, double a[], int m, double am[] )




































{
  int i;
  int j;
  int k;

  k = 0;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      am[k] = a[i];
      k = k + 1;
    }
  }
  return;
}


double *r8vec_stutter_new ( int n, double a[], int m )




































{
  double *am;
  int i;
  int j;
  int k;

  am = new double[m*n];

  k = 0;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      am[k] = a[i];
      k = k + 1;
    }
  }
  return am;
}


double r8vec_sum ( int n, double a[] )































{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }
  return value;
}


double *r8vec_sum_running ( int n, double v[] )




























{
  int i;
  double *s;

  s = new double[n+1];



  s[0] = 0.0;
  for ( i = 1; i < n + 1; i++ )
  {
    s[i] = s[i-1] + v[i-1];
  }

  return s;
}


void r8vec_swap ( int n, double a1[], double a2[] )





























{
  int i;
  double temp;

  for ( i = 0; i < n; i++ )
  {
    temp  = a1[i];
    a1[i] = a2[i];
    a2[i] = temp;
  }

  return;
}


void r8vec_transpose_print ( int n, double a[], string title )








































{
  int i;
  int ihi;
  int ilo;
  int title_length;

  title_length = s_len_trim ( title );

  for ( ilo = 0; ilo < n; ilo = ilo + 5 )
  {
    if ( ilo == 0 )
    {
      cout << title;
    }
    else
    {
      for ( i = 0; i < title_length; i++ )
      {
        cout << " ";
      }
    }
    cout << "  ";
    ihi = i4_min ( ilo + 5, n );
    for ( i = ilo; i < ihi; i++ )
    {
      cout << "  " << setw(12) << a[i];
    }
    cout << "\n";
  }

  return;
}


void r8vec_undex ( int x_num, double x_val[], int x_unique_num, double tol,
  int undx[], int xdnu[] )
































































































{
  int i;
  int *indx;
  int j;



  indx = r8vec_sort_heap_index_a_new ( x_num, x_val );



  i = 0;

  j = 0;
  undx[j] = indx[i];

  xdnu[indx[i]] = j;

  for ( i = 1; i < x_num; i++ )
  {
    if ( tol < fabs ( x_val[indx[i]] - x_val[undx[j]] ) )
    {
      j = j + 1;
      undx[j] = indx[i];
    }
    xdnu[indx[i]] = j;
  }
  delete [] indx;

  return;
}


void r8vec_uniform_01 ( int n, int &seed, double r[] )


































































{
  int i;
  const int i4_huge = 2147483647;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = ( double ) ( seed ) * 4.656612875E-10;
  }

  return;
}


double *r8vec_uniform_01_new ( int n, int &seed )


































































{
  int i;
  const int i4_huge = 2147483647;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = ( double ) ( seed ) * 4.656612875E-10;
  }

  return r;
}


void r8vec_uniform_ab ( int n, double a, double b, int &seed, double x[] )






































































{
  int i;
  const int i4_huge = 2147483647;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    x[i] = a + ( b - a ) * ( double ) ( seed ) * 4.656612875E-10;
  }

  return;
}


double *r8vec_uniform_ab_new ( int n, double a, double b, int &seed )






































































{
  int i;
  const int i4_huge = 2147483647;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_AB_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = a + ( b - a ) * ( double ) ( seed ) * 4.656612875E-10;
  }

  return r;
}


void r8vec_uniform_abvec ( int n, double a[], double b[], int &seed, double x[] )







































































{
  int i;
  const int i4_huge = 2147483647;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_ABVEC - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    x[i] = a[i] + ( b[i] - a[i] ) * ( double ) ( seed ) * 4.656612875E-10;
  }

  return;
}


double *r8vec_uniform_abvec_new ( int n, double a[], double b[], int &seed )








































































{
  int i;
  const int i4_huge = 2147483647;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_ABVEC_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = a[i] + ( b[i] - a[i] ) * ( double ) ( seed ) * 4.656612875E-10;
  }

  return r;
}


double *r8vec_uniform_unit_new ( int m, int &seed )




























{
  double *a;
  int i;
  double norm;



  a = r8vec_normal_01_new ( m, seed );



  norm = 0.0;
  for ( i = 0; i < m; i++ )
  {
    norm = norm + a[i] * a[i];
  }
  norm = sqrt ( norm );



  for ( i = 0; i < m; i++ )
  {
    a[i] = a[i] / norm;
  }

  return a;
}


int r8vec_unique_count ( int n, double a[], double tol )




































{
  int i;
  int j;
  int unique_num;

  unique_num = 0;

  for ( i = 0; i < n; i++ )
  {
    unique_num = unique_num + 1;

    for ( j = 0; j < i; j++ )
    {
      if ( fabs ( a[i] - a[j] ) <= tol )
      {
        unique_num = unique_num - 1;
        break;
      }
    }
  }
  return unique_num;
}


int *r8vec_unique_index ( int n, double a[], double tol )







































{
  int i;
  int j;
  int *unique_index;
  int unique_num;

  unique_index = new int[n];

  for ( i = 0; i < n; i++ )
  {
    unique_index[i] = -1;
  }
  unique_num = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( unique_index[i] == -1 )
    {
      unique_index[i] = unique_num;
      for ( j = i + 1; j < n; j++ )
      {
        if ( fabs ( a[i] - a[j] ) <= tol )
        {
          unique_index[j] = unique_num;
        }
      }
      unique_num = unique_num + 1;
    }
  }
  return unique_index;
}


double r8vec_variance ( int n, double x[] )































{
  int i;
  double mean;
  double variance;

  mean = r8vec_mean ( n, x );

  variance = 0.0;
  for ( i = 0; i < n; i++ )
  {
    variance = variance + ( x[i] - mean ) * ( x[i] - mean );
  }

  if ( 1 < n )
  {
    variance = variance / ( double ) ( n );
  }
  else
  {
    variance = 0.0;
  }

  return variance;
}


double r8vec_variance_circular ( int n, double x[] )
































{
  int i;
  double mean;
  double sum_c;
  double sum_s;
  double value;

  mean = r8vec_mean ( n, x );

  sum_c = 0.0;
  for ( i = 0; i < n; i++ )
  {
    sum_c = sum_c + cos ( x[i] - mean );
  }

  sum_s = 0.0;
  for ( i = 0; i < n; i++ )
  {
    sum_s = sum_s + sin ( x[i] - mean );
  }

  value = sqrt ( sum_c * sum_c + sum_s * sum_s ) / ( double ) n;

  value = 1.0 - value;

  return value;
}


double r8vec_variance_sample ( int n, double x[] )































{
  int i;
  double mean;
  double variance;

  mean = r8vec_mean ( n, x );

  variance = 0.0;
  for ( i = 0; i < n; i++ )
  {
    variance = variance + ( x[i] - mean ) * ( x[i] - mean );
  }

  if ( 1 < n )
  {
    variance = variance / ( double ) ( n - 1 );
  }
  else
  {
    variance = 0.0;
  }

  return variance;
}


void r8vec_variance_update ( int nm1, double mean_nm1, double variance_nm1, 
  double xn, int &n, double &mean_n, double &variance_n )







































{
  if ( nm1 <= 0 )
  {
    n = 1;
    mean_n = xn;
    variance_n = 0.0;
  }
  else
  {
    n = nm1 + 1;
    mean_n = mean_nm1 + ( xn - mean_nm1 ) / ( double ) ( n );
    variance_n = ( variance_nm1 * ( double ) ( nm1 ) 
      + ( xn - mean_nm1 ) * ( xn - mean_n ) ) / ( double ) ( n );
  }

  return;
}


double *r8vec_vector_triple_product ( double v1[3], double v2[3], double v3[3] )


































{
  double *v123;
  double *v23;

  v23 = r8vec_cross_product_3d ( v2, v3 );

  v123 = r8vec_cross_product_3d ( v1, v23 );

  delete [] v23;

  return v123;
}


void r8vec_write ( int n, double r[], string output_file )
































{
  int i;
  ofstream output;

  output.open ( output_file.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8VEC_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }

  for ( i = 0; i < n; i++ )
  {
    output << "  " << setw(16) << r[i] << "\n";
  }

  output.close ( );

  return;
}


void r8vec_zeros ( int n, double a[] )





























{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return;
}


double *r8vec_zeros_new ( int n )





























{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return a;
}


int r8vec2_compare ( int n, double a1[], double a2[], int i, int j )






































{
  int isgn;

  isgn = 0;

  if ( a1[i-1] < a1[j-1] )
  {
    isgn = -1;
  }
  else if ( a1[i-1] == a1[j-1] )
  {
    if ( a2[i-1] < a2[j-1] )
    {
      isgn = -1;
    }
    else if ( a2[i-1] < a2[j-1] )
    {
      isgn = 0;
    }
    else if ( a2[j-1] < a2[i-1] )
    {
      isgn = +1;
    }
  }
  else if ( a1[j-1] < a1[i-1] )
  {
    isgn = +1;
  }

  return isgn;
}


void r8vec2_print ( int n, double a1[], double a2[], string title )
































{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i <= n - 1; i++ )
  {
    cout << setw(6)  << i
         << ": " << setw(14) << a1[i]
         << "  " << setw(14) << a2[i] << "\n";
  }

  return;
}


void r8vec2_print_some ( int n, double x1[], double x2[], int max_print,
  string title )











































{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      cout << setw(6)  << i << ": "
           << setw(14) << x1[i] << "  "
           << setw(14) << x2[i] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print-2; i++ )
    {
      cout << setw(6)  << i     << ": "
           << setw(14) << x1[i] << "  "
           << setw(14) << x2[i] << "\n";
    }
    cout << "......  ..............  ..............\n";
    i = n - 1;
    cout << setw(6)  << i     << ": "
         << setw(14) << x1[i] << "  "
         << setw(14) << x2[i] << "\n";
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      cout << setw(6)  << i     << ": "
           << setw(14) << x1[i] << "  "
           << setw(14) << x2[i] << "\n";
    }
    i = max_print - 1;
    cout << setw(6)  << i     << ": "
         << setw(14) << x1[i] << "  "
         << setw(14) << x2[i] << "...more entries...\n";
  }

  return;
}


void r8vec2_sort_a ( int n, double a1[], double a2[] )

































{
  int i;
  int indx;
  int isgn;
  int j;
  double temp;



  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;



  for ( ; ; )
  {
    sort_heap_external ( n, indx, i, j, isgn );



    if ( 0 < indx )
    {
      temp    = a1[i-1];
      a1[i-1] = a1[j-1];
      a1[j-1] = temp;

      temp    = a2[i-1];
      a2[i-1] = a2[j-1];
      a2[j-1] = temp;
    }



    else if ( indx < 0 )
    {
      isgn = r8vec2_compare ( n, a1, a2, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }
  }

  return;
}


void r8vec2_sort_d ( int n, double a1[], double a2[] )

































{
  int i;
  int indx;
  int isgn;
  int j;
  double temp;



  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;



  for ( ; ; )
  {
    sort_heap_external ( n, indx, i, j, isgn );



    if ( 0 < indx )
    {
      temp    = a1[i-1];
      a1[i-1] = a1[j-1];
      a1[j-1] = temp;

      temp    = a2[i-1];
      a2[i-1] = a2[j-1];
      a2[j-1] = temp;
    }



    else if ( indx < 0 )
    {
      isgn = - r8vec2_compare ( n, a1, a2, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }
  }

  return;
}


int *r8vec2_sort_heap_index_a ( int n, double x[], double y[] )



























































{
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;
  double xval;
  double yval;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = new int[n];

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    indx[0] = indx[0];
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      xval = x[indxt];
      yval = y[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      xval = x[indxt];
      yval = y[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( x[indx[j-1]] < x[indx[j]] ||
          ( x[indx[j-1]] == x[indx[j]] && y[indx[j-1]] < y[indx[j]] ) )
        {
          j = j + 1;
        }
      }

      if ( xval < x[indx[j-1]] ||
         ( xval == x[indx[j-1]] && yval < y[indx[j-1]] ) )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return indx;
}


void r8vec2_sorted_unique ( int n, double a1[], double a2[], int &unique_num )










































{
  int itest;

  unique_num = 0;

  if ( n <= 0 )
  {
    return;
  }

  unique_num = 1;

  for ( itest = 1; itest < n; itest++ )
  {
    if ( a1[itest] != a1[unique_num-1] ||
         a2[itest] != a2[unique_num-1] )
    {
      a1[unique_num] = a1[itest];
      a2[unique_num] = a2[itest];
      unique_num = unique_num + 1;
    }
  }

  return;
}


void r8vec2_sorted_unique_index ( int n, double a1[], double a2[],
  int &unique_num, int indx[] )














































{
  int itest;

  if ( n <= 0 )
  {
    unique_num = 0;
    return;
  }
  i4vec_zeros ( n, indx );

  unique_num = 1;
  indx[0] = 1;

  for ( itest = 2; itest <= n; itest++ )
  {
    if ( a1[itest-2] != a1[itest-1] || a2[itest-2] != a2[itest-1] )
    {
      unique_num = unique_num + 1;
      indx[unique_num-1] = itest;
    }
  }

  return;
}


int r8vec2_sum_max_index ( int n, double a[], double b[] )

































{
  int i;
  double sum_max;
  int sum_max_index;

  if ( n <= 0 )
  {
    sum_max_index = -1;
  }
  else
  {
    sum_max_index = 1;
    sum_max = a[0] + b[0];

    for ( i = 2; i <= n; i++ )
    {
      if ( sum_max < a[i-1] + b[i-1] )
      {
        sum_max = a[i-1] + b[i-1];
        sum_max_index = i;
      }
    }
  }
  return sum_max_index;
}


void r8vec3_print ( int n, double a1[], double a2[], double a3[], string title )




























{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i <= n - 1; i++ )
  {
    cout << setw(4)  << i     << ": "
         << setw(10) << a1[i] << "  "
         << setw(10) << a2[i] << "  "
         << setw(10) << a3[i] << "\n";
  }

  return;
}


double *roots_to_r8poly ( int n, double x[] )



























{
  double *c;
  int i;
  int j;

  c = r8vec_zeros_new ( n + 1 );




  c[n] = 1.0;




  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= n+1-j; i++ )
    {
      c[n-i] = c[n-i] - x[n+1-i-j] * c[n-i+1];
    }
  }
  return c;
}


int s_len_trim ( string s )


























{
  int n;

  n = s.length ( );

  while ( 0 < n )
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}


void sort_heap_external ( int n, int &indx, int &i, int &j, int isgn )





























































{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;



  if ( indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }



  else if ( indx < 0 )
  {
    if ( indx == -2 )
    {
      if ( isgn < 0 )
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      indx = -1;
      i = i_save;
      j = j_save;
      return;
    }

    if ( 0 < isgn )
    {
      indx = 2;
      i = i_save;
      j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      if ( n1 == 1 )
      {
        i_save = 0;
        j_save = 0;
        indx = 0;
      }
      else
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        indx = 1;
      }
      i = i_save;
      j = j_save;
      return;
    }
    k = k - 1;
    k1 = k;
  }



  else if ( indx == 1 )
  {
    k1 = k;
  }

  for ( ; ; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 )
    {
      j_save = k1;
      k1 = i_save;
      indx = -1;
      i = i_save;
      j = j_save;
      return;
    }
    else if ( i_save <= n1 )
    {
      j_save = i_save + 1;
      indx = -2;
      i = i_save;
      j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 )
  {
    i_save = 0;
    j_save = 0;
    indx = 0;
    i = i_save;
    j = j_save;
  }
  else
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    indx = 1;
    i = i_save;
    j = j_save;
  }

  return;
}


void timestamp ( )



























{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
