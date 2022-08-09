# include <iostream>
# include <vector>
# include <numeric>
# include <cmath>
# include <algorithm>
# include <assert.h>
# include <tuple>

const double tol = 1e-4;

// ===================================================================

template <typename T>
void MyPrint(std::vector<T> arr)
{
    for (int i=0; i<arr.size(); i++)
    {
        std::cout << arr[i] << " " ;
    }
    std::cout << std::endl;
}

template <typename T>
std::vector<T> Slice(std::vector<T> arr1, int head, int tail)
{
    std::vector<T> arr2;
    typename std::vector<T>::const_iterator First;
    typename std::vector<T>::const_iterator Second;
    First = arr1.begin() + head;
    Second = arr1.begin() + tail;
    //std::vector<T>::const_iterator First = arr1.begin() + head;
    //std::vector<T>::const_iterator Second = arr1.begin() + tail;
    arr2.assign(First, Second);
    return arr2;
}

template <typename T>
T get_inner_product(std::vector<T> v1, std::vector<T> v2)
{
    std::cout << "v1 and v2:" << std::endl;
    MyPrint(v1);
    MyPrint(v2);
    T dot = std::inner_product(std::begin(v1), std::end(v1), std::begin(v2), 0.);
    std::cout << "dot: " << dot << std::endl;
    return dot;
}

template <typename T>
std::vector<T> scalar_product(std::vector<T> v1, T a)
{
    std::vector<T> v2;
    for (int i=0; i < v1.size(); i++)
    {
        v2.push_back(v1[i] * a);
    }
    return v2;
}

template <typename T>
std::vector<T> vector_add(std::vector<T> v1, std::vector<T> v2)
{
    std::vector<T> result(v1.size(), 0);
    transform(v1.begin(), v1.end(), v2.begin(), result.begin(), std::plus<T>());
    return result;
}


template <typename T>
bool vector_equal(std::vector<T> v1, std::vector<T> v2)
{
    bool tr = 1;
    for (int i = 0; i < v1.size(); i++)
    {
        if (abs(v1[i] - v2[i]) > tol)
        {
            tr = 0;
            break;
        }
    }
    return (tr && v1.size() == v2.size());
}


template <typename T>
bool vector_in(std::vector<std::vector<T>> vv, std::vector<T> v)
{
    bool tr = 0;
    for (int i = 0; i < vv.size(); i++)
    {
        if (vector_equal(vv[i], v))
        {
            tr = 1;
            break;
        }
    }
    return tr;
}



// ===================================================================

//void get_segment(std::vector<std::vector<double>> ineq2, int iline)
//std::tuple<double, double> get_segment(std::vector<std::vector<double>> ineq2, int iline)
std::tuple<std::vector<double>, std::vector<double>> get_segment(std::vector<std::vector<double>> ineq2, int iline)
{
    bool exist = 1;

    auto a0 = scalar_product( Slice(ineq2[iline], 0, 2), ineq2[iline][2] );
    a0 = scalar_product(a0, 1/get_inner_product(Slice(ineq2[iline], 0, 2), Slice(ineq2[iline], 0, 2)));
    //std::cout << "a0:" << std::endl;
    //MyPrint(a0);

    std::vector<double> a1;
    a1.push_back(ineq2[iline][1]);
    a1.push_back(-ineq2[iline][0]);
    //std::cout << "a1:" << std::endl;
    //MyPrint(a1); 

    double xmax = 1.0/1e-6;
    double xmin = -xmax;

    for (int ii = 0; ii < ineq2.size(); ii++)
    {
         std::cout << "ii: " << ii << std::endl;
         if (ii == iline)
         {
             continue;
         } 
         double ff = get_inner_product(Slice(ineq2[ii], 0, 2), a1);
         double dd = get_inner_product(Slice(ineq2[ii], 0, 2), a0);

         if (abs(ff) < tol)
         {
             //std::cout << "ff, dd: " << ff << " " << dd << std::endl;
             //std::cout << "ineq2[ii][2] " << ineq2[ii][2] << std::endl;
             exist = (exist  && dd < ineq2[ii][2] + tol);
             if (exist==0)
             {
                 std::cout << "break" << std::endl;
                 break;
             }
         }
         else if (ff > 0.0)
         {
             xmax = std::min(xmax, (ineq2[ii][2] - dd)/ff);
             std::cout << "xmax: " << xmax << std::endl;
         }
         else 
         {
             xmin = std::max(xmin, (ineq2[ii][2] - dd)/ff);
             std::cout << "xmin: " << xmin << std::endl;
 
         }
         std::cout << std::endl;
    }
    exist = (exist && xmin + tol < xmax);     
    //std::cout << "exist: " << exist << std::endl; 
    
    if (exist)
    {
        assert (xmin > -1.0/0.0 && xmax < 1.0/0.0);
        std::vector<double> vtx1, vtx2;
        vtx1 = vector_add(a0, scalar_product(a1, xmin));
        vtx2 = vector_add(a0, scalar_product(a1, xmax));
        std::tuple<std::vector<double>, std::vector<double>> seg = std::make_tuple(vtx1, vtx2);
        return seg;
    }
    else
    {
        throw "no existing";
    }
  
    //std::tuple<double, double> tu = std::make_tuple(xmin, xmax);
    //return tu;    
}


//void get_polygen(std::vector<std::vector<double>>ineq3, int iplane)
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<int>>, std::vector<int>> get_polygen(std::vector<std::vector<double>>ineq3, int iplane)
{
    //
    // Any point on the selected plane can be written as a0 + x1*a1+x2*a2
    //
    auto a0 = scalar_product(Slice(ineq3[iplane], 0, 3), ineq3[iplane][3]);
    a0 = scalar_product(a0, 1/get_inner_product(Slice(ineq3[iplane], 0, 3), 
                                                Slice(ineq3[iplane], 0, 3)) );    

    double ll, mm, nn;
    ll = -ineq3[iplane][0];
    mm = -ineq3[iplane][1];
    nn = -ineq3[iplane][2];

    std::vector<double> a1, a2;

    if (abs(ll)>tol)
    {
        a1 = {-mm, ll, 0.};
        a2 = {-nn, 0., ll};
    }
    else if (abs(mm)>tol)
    {
        a1 = {mm, -ll, 0.};
        a2 = {0., -nn, mm};
    }
    else if (abs(nn)>tol)
    {
        a1 = {nn, 0., -ll};
        a2 = {0., nn, -mm};
    }

    // 2D inequilities projected on the selected plane
    //   ineq3[i,012]*a1*x1 + ineq3[i,012]*a2*x2 < ineq3[i,3] - ineq3[i,012]*a0
    std::vector<std::vector<double>> ineq2;
    
    for (int ii=0; ii < ineq3.size(); ii++)
    {
        if (ii==iplane) 
        {
            continue;
        }
        double aa1 = get_inner_product(Slice(ineq3[ii], 0, 3), a1);
        double aa2 = get_inner_product(Slice(ineq3[ii], 0, 3), a2);
        double bb = ineq3[ii][3] - get_inner_product(Slice(ineq3[ii], 0, 3), a0);
        // aa1*x1 + aa2*x2 < bb
        if (abs(aa1)>tol || abs(aa2)>tol)
        {
            std::vector<double> tmp = {aa1, aa2, bb};
            ineq2.push_back(tmp);
        }
        else if (0 < bb+tol)
        {
            continue;
        }
        else
        {
            throw "ineq2 wrong !";
        }
    }
 
    // The polygen
    
    std::vector<std::vector<double>> vertex;
    std::vector<std::vector<int>> seg;

    for (int iline=0; iline < ineq2.size(); iline++)
    {
        auto v1_v2 = get_segment(ineq2, iline);
        auto v1_ = std::get<0>(v1_v2);
        auto v2_ = std::get<1>(v1_v2);
 
        //v1 = lattvec(a0+v1_[0]*a1+v1_[1]*a2)
        auto v1 = vector_add( a0, vector_add( scalar_product(a1, v1_[0]), scalar_product(a2, v1_[1]) ) );
        //v2=lattvec(a0+v2_[0]*a1+v2_[1]*a2)
        auto v2 = vector_add( a0, vector_add( scalar_product(a1, v2_[0]), scalar_product(a2, v2_[1]) ) );
          
        int n1;
        if (!vector_in(vertex, v1))
        { 
            vertex.push_back(v1);
            n1 = vertex.size() - 1;
        }
        else
        {
            //n1 = [i for i,v in enumerate(vertex) if v==v1][0];
            n1 = 0;
            while( !vector_equal(vertex[n1], v1) )
            {
                n1++;
            }
        }

        int n2;
        if (!vector_in(vertex, v2))
        {
            vertex.push_back(v2);
            n2 = vertex.size() - 1;
        }
        else
        {
            n2 = 0;
            while( !vector_equal( vertex[n2], v2 ) )
            {
                n2++;
            }
        }
        // add the seg 
        std::vector<int> ee;
        if (n1 <= n2)
        {
            ee = {n1, n2};
        }
        else
        {
            ee = {n2, n1};
        }
        if (!vector_in(seg, ee))
        {
            seg.push_back(ee);
        }
    }  

    // check the connectivity of the polygen
    
    //v2seg=[[] for ii in range(len(vertex)) ]
    std::vector<std::vector<int>> v2seg(vertex.size());
    for (int ii=0; ii<seg.size(); ii++)
    {
        auto ee = seg[ii];
        v2seg[ee[0]].push_back(ii);
        v2seg[ee[1]].push_back(ii);
    }

    bool closed = std::all_of(v2seg.begin(), v2seg.end(), [](std::vector<int> v) 
                                                         { return v.size() == 2; });
 
    closed = (closed && vertex.size() > 0);
    if (!closed)
    {
        throw "polygen not closed !";
    }

    // sort the vertices in the polygen

    std::vector<int> polygen = {0};
    std::vector<bool> passed(seg.size(), 0);  

    for (int ii=0; ii<vertex.size()-1; ii++)
    {
        int ee = v2seg[polygen[-1]][0];  // ee is a seg
        if (passed[ee])
        {
            ee = v2seg[polygen[-1]][1];
        }
        assert (!passed[ee]);

        if ( std::count(polygen.begin(), polygen.end(), seg[ee][0]) 
            && (!std::count(polygen.begin(), polygen.end(), seg[ee][1])) )
        {
            polygen.push_back(seg[ee][1]);
        }
        else if ( std::count(polygen.begin(), polygen.end(), seg[ee][1]) 
            && (!std::count(polygen.begin(), polygen.end(), seg[ee][0])) )
        {
            polygen.push_back(seg[ee][0]);
        }
        else
        {
            throw "point error in polygen !";
        }
        passed[ee] = 1;
    }
    auto tu = std::make_tuple(vertex, seg, polygen);
    return tu;
}



// ===================================================================

void test_print()
{
    std::vector<int> vv = {1,1,1,3};
    MyPrint(vv);

}


void test_slice()
{
    std::vector<int> vv = {1,2,3,4,5,6};
    MyPrint(Slice(vv, 0, 2));
}


void test_get_inner_product()
{
    std::vector<int> vv = {1,2,3};
    std::cout << get_inner_product(vv, vv) << std::endl;
}


void test_vector_add()
{
    std::vector<int> v1 = {1,2,3};
    std::vector<int> v2 = {1,2,3};
    std::vector<int> v3 = vector_add(v1, v2);
    MyPrint(v3);
}


void test_vector_equal()
{
    std::vector<int> v1 = {1,2,3};
    std::vector<int> v2 = {1,2,3,5};
    std::cout << vector_equal(v1, v2) << std::endl; 
}


void test_vector_in()
{
    std::vector<int> v = {1,2,3};
    std::vector<std::vector<int>> vv = {{1,2,3}, {4,5,6}};
    std::cout << vector_in(vv, v) << std::endl; 
}

void test_get_segment()
{
    std::vector<std::vector<double>> ineq2 = {{1, 1, 1}, {-1, 1, 1}, {-1, -1, 1}, {1, -1, 1}};
    auto seg = get_segment(ineq2, 0);
    //std::vector<double> xmin, xmax;
    auto v1 = std::get<0>(seg);
    auto v2= std::get<1>(seg);
    //std::cout << xmin << ", " << xmax << std::endl;
    MyPrint(v1);
    MyPrint(v2);
}


// ====================================================================

int main()
{
    //test_print();
    //test_slice();
    //test_get_inner_product();
    //test_vector_add();
    //test_vector_equal();
    test_vector_in();
    //test_get_segment();
}

