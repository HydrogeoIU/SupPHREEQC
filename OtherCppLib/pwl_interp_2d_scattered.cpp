# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iostream>
# include <iomanip>

using namespace std;

# include "pwl_interp_2d_scattered.hpp"
# include "r8lib.hpp"



int diaedg ( double x0, double y0, double x1, double y1, double x2, double y2,
  double x3, double y3 )













































{
  double ca;
  double cb;
  double dx10;
  double dx12;
  double dx30;
  double dx32;
  double dy10;
  double dy12;
  double dy30;
  double dy32;
  double s;
  double tol;
  double tola;
  double tolb;
  int value;

  tol = 100.0 * r8_epsilon ( );

  dx10 = x1 - x0;
  dy10 = y1 - y0;
  dx12 = x1 - x2;
  dy12 = y1 - y2;
  dx30 = x3 - x0;
  dy30 = y3 - y0;
  dx32 = x3 - x2;
  dy32 = y3 - y2;

  tola = tol * r8_max ( fabs ( dx10 ),
               r8_max ( fabs ( dy10 ),
               r8_max ( fabs ( dx30 ), fabs ( dy30 ) ) ) );

  tolb = tol * r8_max ( fabs ( dx12 ),
               r8_max ( fabs ( dy12 ),
               r8_max ( fabs ( dx32 ), fabs ( dy32 ) ) ) );

  ca = dx10 * dx30 + dy10 * dy30;
  cb = dx12 * dx32 + dy12 * dy32;

  if ( tola < ca && tolb < cb )
  {
    value = -1;
  }
  else if ( ca < -tola && cb < -tolb )
  {
    value = 1;
  }
  else
  {
    tola = r8_max ( tola, tolb );
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb
      + ( dx32 * dy12 - dx12 * dy32 ) * ca;

    if ( tola < s )
    {
      value = -1;
    }
    else if ( s < -tola )
    {
      value = 1;
    }
    else
    {
      value = 0;
    }

  }

  return value;
}


void i4mat_transpose_print ( int m, int n, int a[], string title )

































{
  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}


void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, string title )






































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



  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    cout << "\n";





    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i - 1 << "  ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";



    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {



      cout << setw(5) << j - 1 << ":";
      for ( i = i2lo; i <= i2hi; i++ )
      {
        cout << setw(6) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}


void i4vec_heap_d ( int n, int a[] )




















































{
  int i;
  int ifree;
  int key;
  int m;



  for ( i = ( n / 2 ) - 1; 0 <= i; i-- )
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


int i4vec_min ( int n, int a[] )
































{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] < value )
    {
      value = a[i];
    }
  }
  return value;
}


void i4vec_sort_heap_a ( int n, int a[] )






































{
  int n1;
  int temp;

  if ( n <= 1 )
  {
    return;
  }



  i4vec_heap_d ( n, a );






  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;



  for ( n1 = n-1; 2 <= n1; n1-- )
  {



    i4vec_heap_d ( n1, a );



    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;
  }
  return;
}


int i4vec_sorted_unique ( int n, int a[] )
































{
  int i;
  int unique_num;

  unique_num = 0;

  if ( n <= 0 )
  {
    return unique_num;
  }

  unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] != a[unique_num-1] )
    {
      unique_num = unique_num + 1;
      a[unique_num-1] = a[i];
    }
  }

  return unique_num;
}


int lrline ( double xu, double yu, double xv1, double yv1, double xv2,
  double yv2, double dv )














































{
  double dx;
  double dxu;
  double dy;
  double dyu;
  double t;
  double tol = 0.0000001;
  double tolabs;
  int value;

  dx = xv2 - xv1;
  dy = yv2 - yv1;
  dxu = xu - xv1;
  dyu = yu - yv1;

  tolabs = tol * r8_max ( fabs ( dx ),
                 r8_max ( fabs ( dy ),
                 r8_max ( fabs ( dxu ),
                 r8_max ( fabs ( dyu ), fabs ( dv ) ) ) ) );

  t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy );

  if ( tolabs < t )
  {
    value = 1;
  }
  else if ( -tolabs <= t )
  {
    value = 0;
  }
  else if ( t < -tolabs )
  {
    value = -1;
  }

  return value;
}


int perm_check2 ( int n, int p[], int base )





































{
  int error;
  int ifind;
  int iseek;

  error = 0;

  for ( iseek = base; iseek < base + n; iseek++ )
  {
    error = 1;

    for ( ifind = 1; ifind <= n; ifind++ )
    {
      if ( p[ifind-1] == iseek )
      {
        error = 0;
        break;
      }
    }

    if ( error )
    {
      cerr << "\n";
      cerr << "PERM_CHECK2 - Fatal error!\n";
      cerr << "  Could not find occurrence of value " << iseek << "\n";
      return 1;
    }
  }

  return 0;
}


void perm_inverse ( int n, int p[] )


























{
  int base;
  int error;
  int i;
  int i0;
  int i1;
  int i2;
  int is;
  int p_min;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "PERM_INVERSE - Fatal error!\n";
    cerr << "  Input value of N = " << n << "\n";
    exit ( 1 );
  }



  p_min = i4vec_min ( n, p );
  base = 1;

  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - p_min + base;
  }



  error = perm_check2 ( n, p, base );

  if ( error )
  {
    cerr << "\n";
    cerr << "PERM_INVERSE - Fatal error!\n";
    cerr << "  The input array does not represent\n";
    cerr << "  a proper permutation.\n";
    exit ( 1 );
  }

  is = 1;

  for ( i = 1; i <= n; i++ )
  {
    i1 = p[i-1];

    while ( i < i1 )
    {
      i2 = p[i1-1];
      p[i1-1] = -i2;
      i1 = i2;
    }

    is = - i4_sign ( p[i-1] );
    p[i-1] = abs ( p[i-1] ) * i4_sign ( is );

  }

  for ( i = 1; i <= n; i++ )
  {
    i1 = -p[i-1];

    if ( 0 <= i1 )
    {
      i0 = i;

      for ( ; ; )
      {
        i2 = p[i1-1];
        p[i1-1] = i0;

        if ( i2 < 0 )
        {
          break;
        }

        i0 = i1;
        i1 = i2;
      }
    }
  }



  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + p_min - base;
  }

  return;
}


double *pwl_interp_2d_scattered_value ( int nd, double xyd[], double zd[], 
  int t_num, int t[], int t_neighbor[], int ni, double xyi[] )







































{
  double alpha;
  double beta;
  double gamma;
  int edge;
  int i;
  int j;
  int step_num;
  double *zi;

  zi = new double[ni];

  for ( i = 0; i < ni; i++ )
  {
    triangulation_search_delaunay ( nd, xyd, 3, t_num, t, t_neighbor, 
      xyi+2*i, j, alpha, beta, gamma, edge, step_num );

    if ( j == -1 )
    {
      zi[i] = -1.0;
    }

    zi[i] = alpha * zd[t[0+j*3]] 
          + beta  * zd[t[1+j*3]] 
          + gamma * zd[t[2+j*3]];
  }
  return zi;
}


int r8tris2 ( int node_num, double node_xy[], int &triangle_num,
  int triangle_node[], int triangle_neighbor[] )
























































{
  int base;
  double cmax;
  int e;
  int error;
  int i;
  int *indx;
  int j;
  int k;
  int l;
  int ledg;
  int lr;
  int ltri;
  int m;
  int m1;
  int m2;
  int n;
  int redg;
  int rtri;
  int *stack;
  int t;
  double tol;
  int top;

  stack = new int[node_num];

  tol = 100.0 * r8_epsilon ( );



  base = 0;

  indx = r82row_sort_heap_index_a ( node_num, node_xy );

  r82row_permute ( node_num, indx, node_xy );



  m1 = 1;

  for ( i = 2; i <= node_num; i++ )
  {
    m = m1;
    m1 = i;

    k = -1;

    for ( j = 0; j <= 1; j++ )
    {
      cmax = r8_max ( fabs ( node_xy[2*(m-1)+j] ),
                     fabs ( node_xy[2*(m1-1)+j] ) );

      if ( tol * ( cmax + 1.0 )
           < fabs ( node_xy[2*(m-1)+j] - node_xy[2*(m1-1)+j] ) )
      {
        k = j;
        break;
      }

    }

    if ( k == -1 )
    {
      cerr << "\n";
      cerr << "R8TRIS2 - Fatal error!\n";
      cerr << "  Fails for point number I = " << i << "\n";
      cerr << "  M =  " << m  << "\n";
      cerr << "  M1 = " << m1 << "\n";
      cerr << "  X,Y(M)  = " << node_xy[2*(m-1)+0] << "  "
                             << node_xy[2*(m-1)+1] << "\n";
      cerr << "  X,Y(M1) = " << node_xy[2*(m1-1)+0] << "  "
                             << node_xy[2*(m1-1)+1] << "\n";
      exit ( 1 );
    }

  }




  m1 = 1;
  m2 = 2;
  j = 3;

  for ( ; ; )
  {
    if ( node_num < j )
    {
      cerr << "\n";
      cerr << "R8TRIS2 - Fatal error!\n";
      delete [] stack;
      return 225;
    }

    m = j;

    lr = lrline ( node_xy[2*(m-1)+0], node_xy[2*(m-1)+1],
      node_xy[2*(m1-1)+0], node_xy[2*(m1-1)+1],
      node_xy[2*(m2-1)+0], node_xy[2*(m2-1)+1], 0.0 );

    if ( lr != 0 )
    {
      break;
    }

    j = j + 1;

  }




  triangle_num = j - 2;

  if ( lr == -1 )
  {
    triangle_node[3*0+0] = m1;
    triangle_node[3*0+1] = m2;
    triangle_node[3*0+2] = m;
    triangle_neighbor[3*0+2] = -3;

    for ( i = 2; i <= triangle_num; i++ )
    {
      m1 = m2;
      m2 = i+1;
      triangle_node[3*(i-1)+0] = m1;
      triangle_node[3*(i-1)+1] = m2;
      triangle_node[3*(i-1)+2] = m;
      triangle_neighbor[3*(i-1)+0] = -3 * i;
      triangle_neighbor[3*(i-1)+1] = i;
      triangle_neighbor[3*(i-1)+2] = i - 1;

    }

    triangle_neighbor[3*(triangle_num-1)+0] = -3 * triangle_num - 1;
    triangle_neighbor[3*(triangle_num-1)+1] = -5;
    ledg = 2;
    ltri = triangle_num;
  }
  else
  {
    triangle_node[3*0+0] = m2;
    triangle_node[3*0+1] = m1;
    triangle_node[3*0+2] = m;
    triangle_neighbor[3*0+0] = -4;

    for ( i = 2; i <= triangle_num; i++ )
    {
      m1 = m2;
      m2 = i+1;
      triangle_node[3*(i-1)+0] = m2;
      triangle_node[3*(i-1)+1] = m1;
      triangle_node[3*(i-1)+2] = m;
      triangle_neighbor[3*(i-2)+2] = i;
      triangle_neighbor[3*(i-1)+0] = -3 * i - 3;
      triangle_neighbor[3*(i-1)+1] = i - 1;
    }

    triangle_neighbor[3*(triangle_num-1)+2] = -3 * triangle_num;
    triangle_neighbor[3*0+1] = -3 * triangle_num - 2;
    ledg = 2;
    ltri = 1;

  }





  top = 0;

  for ( i = j+1; i <= node_num; i++ )
  {
    m = i;
    m1 = triangle_node[3*(ltri-1)+ledg-1];

    if ( ledg <= 2 )
    {
      m2 = triangle_node[3*(ltri-1)+ledg];
    }
    else
    {
      m2 = triangle_node[3*(ltri-1)+0];
    }

    lr = lrline ( node_xy[2*(m-1)+0], node_xy[2*(m-1)+1],
      node_xy[2*(m1-1)+0], node_xy[2*(m1-1)+1],
      node_xy[2*(m2-1)+0], node_xy[2*(m2-1)+1], 0.0 );

    if ( 0 < lr )
    {
      rtri = ltri;
      redg = ledg;
      ltri = 0;
    }
    else
    {
      l = -triangle_neighbor[3*(ltri-1)+ledg-1];
      rtri = l / 3;
      redg = (l % 3) + 1;
    }

    vbedg ( node_xy[2*(m-1)+0], node_xy[2*(m-1)+1], node_num,
      node_xy, triangle_num, triangle_node, triangle_neighbor,
      ltri, ledg, rtri, redg );

    n = triangle_num + 1;
    l = -triangle_neighbor[3*(ltri-1)+ledg-1];

    for ( ; ; )
    {
      t = l / 3;
      e = ( l % 3 ) + 1;
      l = -triangle_neighbor[3*(t-1)+e-1];
      m2 = triangle_node[3*(t-1)+e-1];

      if ( e <= 2 )
      {
        m1 = triangle_node[3*(t-1)+e];
      }
      else
      {
        m1 = triangle_node[3*(t-1)+0];
      }

      triangle_num = triangle_num + 1;
      triangle_neighbor[3*(t-1)+e-1] = triangle_num;
      triangle_node[3*(triangle_num-1)+0] = m1;
      triangle_node[3*(triangle_num-1)+1] = m2;
      triangle_node[3*(triangle_num-1)+2] = m;
      triangle_neighbor[3*(triangle_num-1)+0] = t;
      triangle_neighbor[3*(triangle_num-1)+1] = triangle_num - 1;
      triangle_neighbor[3*(triangle_num-1)+2] = triangle_num + 1;
      top = top + 1;

      if ( node_num < top )
      {
        cerr << "\n";
        cerr << "R8TRIS2 - Fatal error!\n";
        cerr << "  Stack overflow.\n";
        delete [] stack;
        return 8;
      }

      stack[top-1] = triangle_num;

      if ( t == rtri && e == redg )
      {
        break;
      }

    }

    triangle_neighbor[3*(ltri-1)+ledg-1] = -3 * n - 1;
    triangle_neighbor[3*(n-1)+1] = -3 * triangle_num - 2;
    triangle_neighbor[3*(triangle_num-1)+2] = -l;
    ltri = n;
    ledg = 2;

    error = swapec ( m, top, ltri, ledg, node_num, node_xy, triangle_num,
      triangle_node, triangle_neighbor, stack );

    if ( error != 0 )
    {
      cerr << "\n";
      cerr << "R8TRIS2 - Fatal error!\n";
      cerr << "  Error return from SWAPEC.\n";
      delete [] stack;
      return error;
    }

  }



  for ( i = 0; i < 3; i++ )
  {
    for ( j = 0; j < triangle_num; j++ )
    {
      triangle_node[i+j*3] = indx [ triangle_node[i+j*3] - 1 ];
    }
  }

  perm_inverse ( node_num, indx );

  r82row_permute ( node_num, indx, node_xy );

  delete [] indx;
  delete [] stack;

  return 0;
}


int swapec ( int i, int &top, int &btri, int &bedg, int node_num,
  double node_xy[], int triangle_num, int triangle_node[],
  int triangle_neighbor[], int stack[] )



































































{
  int a;
  int b;
  int c;
  int e;
  int ee;
  int em1;
  int ep1;
  int f;
  int fm1;
  int fp1;
  int l;
  int r;
  int s;
  int swap;
  int t;
  int tt;
  int u;
  double x;
  double y;




  x = node_xy[2*(i-1)+0];
  y = node_xy[2*(i-1)+1];

  for ( ; ; )
  {
    if ( top <= 0 )
    {
      break;
    }

    t = stack[top-1];
    top = top - 1;

    if ( triangle_node[3*(t-1)+0] == i )
    {
      e = 2;
      b = triangle_node[3*(t-1)+2];
    }
    else if ( triangle_node[3*(t-1)+1] == i )
    {
      e = 3;
      b = triangle_node[3*(t-1)+0];
    }
    else
    {
      e = 1;
      b = triangle_node[3*(t-1)+1];
    }

    a = triangle_node[3*(t-1)+e-1];
    u = triangle_neighbor[3*(t-1)+e-1];

    if ( triangle_neighbor[3*(u-1)+0] == t )
    {
      f = 1;
      c = triangle_node[3*(u-1)+2];
    }
    else if ( triangle_neighbor[3*(u-1)+1] == t )
    {
      f = 2;
      c = triangle_node[3*(u-1)+0];
    }
    else
    {
      f = 3;
      c = triangle_node[3*(u-1)+1];
    }

    swap = diaedg ( x, y,
      node_xy[2*(a-1)+0], node_xy[2*(a-1)+1],
      node_xy[2*(c-1)+0], node_xy[2*(c-1)+1],
      node_xy[2*(b-1)+0], node_xy[2*(b-1)+1] );

    if ( swap == 1 )
    {
      em1 = i4_wrap ( e - 1, 1, 3 );
      ep1 = i4_wrap ( e + 1, 1, 3 );
      fm1 = i4_wrap ( f - 1, 1, 3 );
      fp1 = i4_wrap ( f + 1, 1, 3 );

      triangle_node[3*(t-1)+ep1-1] = c;
      triangle_node[3*(u-1)+fp1-1] = i;
      r = triangle_neighbor[3*(t-1)+ep1-1];
      s = triangle_neighbor[3*(u-1)+fp1-1];
      triangle_neighbor[3*(t-1)+ep1-1] = u;
      triangle_neighbor[3*(u-1)+fp1-1] = t;
      triangle_neighbor[3*(t-1)+e-1] = s;
      triangle_neighbor[3*(u-1)+f-1] = r;

      if ( 0 < triangle_neighbor[3*(u-1)+fm1-1] )
      {
        top = top + 1;
        stack[top-1] = u;
      }

      if ( 0 < s )
      {
        if ( triangle_neighbor[3*(s-1)+0] == u )
        {
          triangle_neighbor[3*(s-1)+0] = t;
        }
        else if ( triangle_neighbor[3*(s-1)+1] == u )
        {
          triangle_neighbor[3*(s-1)+1] = t;
        }
        else
        {
          triangle_neighbor[3*(s-1)+2] = t;
        }

        top = top + 1;

        if ( node_num < top )
        {
          return 8;
        }

        stack[top-1] = t;
      }
      else
      {
        if ( u == btri && fp1 == bedg )
        {
          btri = t;
          bedg = e;
        }

        l = - ( 3 * t + e - 1 );
        tt = t;
        ee = em1;

        while ( 0 < triangle_neighbor[3*(tt-1)+ee-1] )
        {
          tt = triangle_neighbor[3*(tt-1)+ee-1];

          if ( triangle_node[3*(tt-1)+0] == a )
          {
            ee = 3;
          }
          else if ( triangle_node[3*(tt-1)+1] == a )
          {
            ee = 1;
          }
          else
          {
            ee = 2;
          }

        }

        triangle_neighbor[3*(tt-1)+ee-1] = l;

      }

      if ( 0 < r )
      {
        if ( triangle_neighbor[3*(r-1)+0] == t )
        {
          triangle_neighbor[3*(r-1)+0] = u;
        }
        else if ( triangle_neighbor[3*(r-1)+1] == t )
        {
          triangle_neighbor[3*(r-1)+1] = u;
        }
        else
        {
          triangle_neighbor[3*(r-1)+2] = u;
        }
      }
      else
      {
        if ( t == btri && ep1 == bedg )
        {
          btri = u;
          bedg = f;
        }

        l = - ( 3 * u + f - 1 );
        tt = u;
        ee = fm1;

        while ( 0 < triangle_neighbor[3*(tt-1)+ee-1] )
        {
          tt = triangle_neighbor[3*(tt-1)+ee-1];

          if ( triangle_node[3*(tt-1)+0] == b )
          {
            ee = 3;
          }
          else if ( triangle_node[3*(tt-1)+1] == b )
          {
            ee = 1;
          }
          else
          {
            ee = 2;
          }

        }

        triangle_neighbor[3*(tt-1)+ee-1] = l;
      }
    }
  }

  return 0;
}


void triangulation_order3_print ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] )
















































{
# define DIM_NUM 2

  int boundary_num;
  int i;
  int j;
  int k;
  int n1;
  int n2;
  int s;
  int s1;
  int s2;
  bool skip;
  int t;
  int *vertex_list;
  int vertex_num;

  cout << "\n";
  cout << "TRIANGULATION_ORDER3_PRINT\n";
  cout << "  Information defining a triangulation.\n";
  cout << "\n";
  cout << "  The number of nodes is " << node_num << "\n";

  r8mat_transpose_print ( DIM_NUM, node_num, node_xy, "  Node coordinates" );

  cout << "\n";
  cout << "  The number of triangles is " << triangle_num << "\n";
  cout << "\n";
  cout << "  Sets of three nodes are used as vertices of\n";
  cout << "  the triangles.  For each triangle, the nodes\n";
  cout << "  are listed in counterclockwise order.\n";

  i4mat_transpose_print ( 3, triangle_num, triangle_node, "  Triangle nodes" );

  cout << "\n";
  cout << "  On each side of a given triangle, there is either\n";
  cout << "  another triangle, or a piece of the convex hull.\n";
  cout << "  For each triangle, we list the indices of the three\n";
  cout << "  neighbors, or (if negative) the codes of the\n";
  cout << "  segments of the convex hull.\n";

  i4mat_transpose_print ( 3, triangle_num, triangle_neighbor,
    "  Triangle neighbors" );



  vertex_list = new int[3*triangle_num];

  k = 0;
  for ( t = 0; t < triangle_num; t++ )
  {
    for ( s = 0; s < 3; s++ )
    {
      vertex_list[k] = triangle_node[s+t*3];
      k = k + 1;
    }
  }

  i4vec_sort_heap_a ( 3*triangle_num, vertex_list );

  vertex_num = i4vec_sorted_unique ( 3*triangle_num, vertex_list );

  delete [] vertex_list;



  boundary_num = 2 * vertex_num - triangle_num - 2;

  cout << "\n";
  cout << "  The number of boundary points is " << boundary_num << "\n";
  cout << "\n";
  cout << "  The segments that make up the convex hull can be\n";
  cout << "  determined from the negative entries of the triangle\n";
  cout << "  neighbor list.\n";
  cout << "\n";
  cout << "     #   Tri  Side    N1    N2\n";
  cout << "\n";

  skip = false;

  k = 0;

  for ( i = 0; i < triangle_num; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      if ( triangle_neighbor[j+i*3] < 0 )
      {
        s = -triangle_neighbor[j+i*3];
        t = s / 3;

        if ( t < 1 || triangle_num < t )
        {
          cout << "\n";
          cout << "  Sorry, this data does not use the R8TRIS2\n";
          cout << "  convention for convex hull segments.\n";
          skip = true;
          break;
        }

        s1 = ( s % 3 ) + 1;
        s2 = i4_wrap ( s1+1, 1, 3 );
        k = k + 1;
        n1 = triangle_node[s1-1+(t-1)*3];
        n2 = triangle_node[s2-1+(t-1)*3];
        cout                  << "  "
             << setw(4) << k  << "  "
             << setw(4) << t  << "  "
             << setw(4) << s1 << "  "
             << setw(4) << n1 << "  "
             << setw(4) << n2 << "\n";
      }
    }

    if ( skip )
    {
      break;
    }
  }

  return;
# undef DIM_NUM
}


void triangulation_search_delaunay ( int node_num, double node_xy[],
  int triangle_order, int triangle_num, int triangle_node[],
  int triangle_neighbor[], double p[2], int &triangle_index, 
  double &alpha, double &beta, double &gamma, int &edge,
  int &step_num )




















































































{
  int a;
  int b;
  int c;
  double det;
  double dxp;
  double dxa;
  double dxb;
  double dyp;
  double dya;
  double dyb;
  static int triangle_index_save = -1;

  step_num = - 1;
  edge = 0;

  if ( triangle_index_save < 0 || triangle_num <= triangle_index_save )
  {
    triangle_index = ( triangle_num + 1 ) / 2;
  }
  else
  {
    triangle_index = triangle_index_save;
  }

  for ( ; ; )
  {
    step_num = step_num + 1;

    if ( triangle_num < step_num )
    {
      cerr << "\n";
      cerr << "TRIANGULATION_SEARCH_DELAUNAY - Fatal error!\n";
      cerr << "  The algorithm seems to be cycling.\n";
      cerr << "  Current triangle is " << triangle_index << "\n";
      exit ( 1 );
    }



    a = triangle_node[0+triangle_index*triangle_order];
    b = triangle_node[1+triangle_index*triangle_order];
    c = triangle_node[2+triangle_index*triangle_order];




    dxa = node_xy[0+a*2] - node_xy[0+c*2];
    dya = node_xy[1+a*2] - node_xy[1+c*2];

    dxb = node_xy[0+b*2] - node_xy[0+c*2];
    dyb = node_xy[1+b*2] - node_xy[1+c*2];

    dxp = p[0]           - node_xy[0+c*2];
    dyp = p[1]           - node_xy[1+c*2];

    det = dxa * dyb - dya * dxb;




    alpha = ( dxp * dyb - dyp * dxb ) / det;
    beta =  ( dxa * dyp - dya * dxp ) / det;
    gamma = 1.0 - alpha - beta;




    if ( 0.0 <= alpha &&
         0.0 <= beta  &&
         0.0 <= gamma )
    {
      break;
    }










    if ( alpha < 0.0 && 0 <= triangle_neighbor[1+triangle_index*3] )
    {
      triangle_index = triangle_neighbor[1+triangle_index*3];
      continue;
    }
    else if ( beta < 0.0 && 0 <= triangle_neighbor[2+triangle_index*3] )
    {
      triangle_index = triangle_neighbor[2+triangle_index*3];
      continue;
    }
    else if ( gamma < 0.0 && 0 <= triangle_neighbor[0+triangle_index*3] )
    {
      triangle_index = triangle_neighbor[0+triangle_index*3];
      continue;
    }






    if ( alpha < 0.0 )
    {
      edge = -2;
      break;
    }
    else if ( beta < 0.0 )
    {
      edge = -3;
      break;
    }
    else if ( gamma < 0.0 )
    {
      edge = -1;
      break;
    }
    else
    {
      cerr << "\n";
      cerr << "TRIANGULATION_ORDER3_SEARCH - Fatal error!\n";
      cerr << "  The algorithm seems to have reached a dead end\n";
      cerr << "  after " << step_num << " steps.\n";
      exit ( 1 );
    }
  }
  triangle_index_save = triangle_index;

  return;
}


void vbedg ( double x, double y, int node_num, double node_xy[],
  int triangle_num, int triangle_node[], int triangle_neighbor[], int &ltri,
  int &ledg, int &rtri, int &redg )

































































{
  int a;
  double ax;
  double ay;
  int b;
  double bx;
  double by;
  bool done;
  int e;
  int l;
  int lr;
  int t;




  if ( ltri == 0 )
  {
    done = false;
    ltri = rtri;
    ledg = redg;
  }
  else
  {
    done = true;
  }

  for ( ; ; )
  {
    l = -triangle_neighbor[3*(rtri-1)+redg-1];
    t = l / 3;
    e = 1 + l % 3;
    a = triangle_node[3*(t-1)+e-1];

    if ( e <= 2 )
    {
      b = triangle_node[3*(t-1)+e];
    }
    else
    {
      b = triangle_node[3*(t-1)+0];
    }

    ax = node_xy[2*(a-1)+0];
    ay = node_xy[2*(a-1)+1];

    bx = node_xy[2*(b-1)+0];
    by = node_xy[2*(b-1)+1];

    lr = lrline ( x, y, ax, ay, bx, by, 0.0 );

    if ( lr <= 0 )
    {
      break;
    }

    rtri = t;
    redg = e;

  }

  if ( done )
  {
    return;
  }

  t = ltri;
  e = ledg;

  for ( ; ; )
  {
    b = triangle_node[3*(t-1)+e-1];
    e = i4_wrap ( e-1, 1, 3 );

    while ( 0 < triangle_neighbor[3*(t-1)+e-1] )
    {
      t = triangle_neighbor[3*(t-1)+e-1];

      if ( triangle_node[3*(t-1)+0] == b )
      {
        e = 3;
      }
      else if ( triangle_node[3*(t-1)+1] == b )
      {
        e = 1;
      }
      else
      {
        e = 2;
      }

    }

    a = triangle_node[3*(t-1)+e-1];
    ax = node_xy[2*(a-1)+0];
    ay = node_xy[2*(a-1)+1];

    bx = node_xy[2*(b-1)+0];
    by = node_xy[2*(b-1)+1];

    lr = lrline ( x, y, ax, ay, bx, by, 0.0 );

    if ( lr <= 0 )
    {
      break;
    }

  }

  ltri = t;
  ledg = e;

  return;
}
