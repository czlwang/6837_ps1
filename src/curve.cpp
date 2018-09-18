#include "curve.h"
#include "vertexrecorder.h"
using namespace std;

const float c_pi = 3.14159265358979323846f;

namespace
{
// Approximately equal to.  We don't want to use == because of
// precision issues with floating point.
inline bool approx(const Vector3f& lhs, const Vector3f& rhs)
{
	const float eps = 1e-8f;
	return (lhs - rhs).absSquared() < eps;
}

static bool checkXY(const vector< Vector3f >& P)
{//question: does something still live in the xy plane after being transformed. and vice versa. 
    std::cout << "checking for xy plane " << std::endl;
    for (Vector3f point : P)
    {
        if(point[2] != 0)
        {
            return false;
        }
    }
    return true;
}
    
}

//globals
//Bezier matrix
Matrix4f B(1.0, -3.0,  3.0, -1.0,
           0.0,  3.0, -6.0,  3.0,
           0.0,  0.0,  3.0, -3.0,
           0.0,  0.0,  0.0,  1.0);
//Bspline
Matrix4f BS6(1.0, -3.0,  3.0, -1.0,
             4.0,  0.0, -6.0,  3.0,
             1.0,  3.0,  3.0, -3.0,
             0.0,  0.0,  0.0,  1.0);
Matrix4f BS = (1.0/6.)*BS6;


Curve evalBezier(const vector< Vector3f >& P, unsigned steps)
{
	// Check
	if (P.size() < 4 || P.size() % 3 != 1)
	{
		cerr << "evalBezier must be called with 3n+1 control points." << endl;
		exit(0);
	}
    
    Curve* curve = new std::vector< CurvePoint >();
    
    bool liesInXY = checkXY(P);

	// TODO:
	// You should implement this function so that it returns a Curve
	// (e.g., a vector< CurvePoint >).  The variable "steps" tells you
	// the number of points to generate on each piece of the spline.
	// At least, that's how the sample solution is implemented and how
	// the SWP files are written.  But you are free to interpret this
	// variable however you want, so long as you can control the
	// "resolution" of the discretized spline curve with it.

	// Make sure that this function computes all the appropriate
	// Vector3fs for each CurvePoint: V,T,N,B.
	// [NBT] should be unit and orthogonal.

	// Also note that you may assume that all Bezier curves that you
	// receive have G1 continuity.  Otherwise, the TNB will not be
	// be defined at points where this does not hold.

	cerr << "\t>>> evalBezier has been called with the following input:" << endl;

//	cerr << "\t>>> Control points (type vector< Vector3f >): " << endl;
//	for (int i = 0; i < (int)P.size(); ++i)
//	{
//        cerr << "\t>>> " << P[i] << endl;
//        cerr << "\t>>> " << P[i][0] << endl;
//        cerr << "\t>>> " << P[i][1] << endl;
//        cerr << "\t>>> " << P[i][2] << endl;
//	}

    //velocity matrx
    Matrix4f B_vel(-3.0,   6.0, -3.0, 0.0,
                    3.0, -12.0,  9.0, 0.0,
                    0.0,   6.0, -9.0, 0.0,
                    0.0,   0.0,  3.0, 0.0);
    
    int pieces = (int)(P.size()-1)/3;
    for (int i = 0; i < pieces; i++)
    {
//        std::cout << "iteration " << i << std::endl;
//        std::cout << "P size " << P.size() << std::endl;

        Vector4f p1(P[3*i], 0.0);
        Vector4f p2(P[3*i+1], 0.0);
        Vector4f p3(P[3*i+2], 0.0);
        Vector4f p4(P[3*i+3], 0.0);
        Matrix4f G(p1, p2, p3, p4, true);//geometry matrix - set of points
        
        Vector3f binormal(0,0,1);
        
        for (int j=0; j <= (int)steps; j++)
        {
            float t = ((float)j)/steps;
            Vector4f Monomial(1, t, pow(t, 2), pow(t, 3));
            Vector3f V = (G*B*Monomial).xyz();//vertex
            Vector3f T = (G*B_vel*Monomial).xyz();//tangent
            T.normalize();
//            std::cout << "T " << T[0] << " " << T[1] << " " << T[2] << std::endl;
            
            Vector3f N;
            N = Vector3f::cross(binormal, T);
            N.normalize();
//            std::cout << "N " << N[0] << " " << N[1] << " " << N[2] << std::endl;
            
            binormal = Vector3f::cross(T, N);
            binormal.normalize();
            
            if(liesInXY)
            {
                std::cout << "lies in xy" << std::endl;
                binormal = Vector3f(0,0,1);
                N = Vector3f::cross(binormal, T);
                N.normalize();
            }
//            std::cout << "binormal " << binormal[0] << " " << binormal[1] << " " << binormal[2] << std::endl;
            
            CurvePoint* point = new CurvePoint();
            point->V = V;
            point->T = T;
            point->B = binormal;
            point->N = N;
            curve->push_back(*point);
        }

    }
    
    //loop through and get the points 4 at a time.
    //get a 3x4 matrix of control points - G
    //spline matrix - 4x4 - B
    //multiply by T for each t in time steps

	cerr << "\t>>> Steps (type steps): " << steps << endl;
	cerr << "\t>>> Returning empty curve." << endl;

	// Right now this will just return this empty curve.
	return *curve;
}

Curve evalBspline(const vector< Vector3f >& P, unsigned steps)
{
	// Check
	if (P.size() < 4)
	{
		cerr << "evalBspline must be called with 4 or more control points." << endl;
		exit(0);
	}
    
    Curve* curve = new std::vector< CurvePoint >();

	// TODO:
	// It is suggested that you implement this function by changing
	// basis from B-spline to Bezier.  That way, you can just call
	// your evalBezier function.

	cerr << "\t>>> evalBSpline has been called with the following input:" << endl;

//	cerr << "\t>>> Control points (type vector< Vector3f >): " << endl;
//	for (int i = 0; i < (int)P.size(); ++i)
//	{
//		cerr << "\t>>> " << P[i] << endl;
//	}

    int pieces = (int)P.size()-3;
    for (int i = 0; i < pieces; i++)
    {
//        std::cout << "iteration " << i << std::endl;
//        std::cout << "P size " << P.size() << std::endl;
        
        Vector4f p1(P[i], 0.0);
        Vector4f p2(P[i+1], 0.0);
        Vector4f p3(P[i+2], 0.0);
        Vector4f p4(P[i+3], 0.0);
        Matrix4f G(p1, p2, p3, p4, true);//geometry matrix - set of points
        G = G*BS*B.inverse();
        vector< Vector3f > control_points;
        
        for(int k = 0; k < 4; k++)
        {
            control_points.push_back(G.getCol(k).xyz());
        }
        
        Curve bezier_curve = evalBezier(control_points, steps);
//        std::cout << "bezier curve size " << curve->size() << std::endl;
        curve->insert(curve->end(), bezier_curve.begin(), bezier_curve.end());
    }

	cerr << "\t>>> Steps (type steps): " << steps << endl;
	cerr << "\t>>> Returning empty curve." << endl;

	// Return an empty curve right now.
//    std::cout << "curve size " << curve->size() << std::endl;
	return *curve;
}

Curve evalCircle(float radius, unsigned steps)
{
	// This is a sample function on how to properly initialize a Curve
	// (which is a vector< CurvePoint >).

	// Preallocate a curve with steps+1 CurvePoints
	Curve R(steps + 1);

	// Fill it in counterclockwise
	for (unsigned i = 0; i <= steps; ++i)
	{
		// step from 0 to 2pi
		float t = 2.0f * c_pi * float(i) / steps;

		// Initialize position
		// We're pivoting counterclockwise around the y-axis
		R[i].V = radius * Vector3f(cos(t), sin(t), 0);

		// Tangent vector is first derivative
		R[i].T = Vector3f(-sin(t), cos(t), 0);

		// Normal vector is second derivative
		R[i].N = Vector3f(-cos(t), -sin(t), 0);

		// Finally, binormal is facing up.
		R[i].B = Vector3f(0, 0, 1);
	}

	return R;
}

void recordCurve(const Curve& curve, VertexRecorder* recorder)
{
	const Vector3f WHITE(1, 1, 1);
	for (int i = 0; i < (int)curve.size() - 1; ++i)
	{
		recorder->record_poscolor(curve[i].V, WHITE);
		recorder->record_poscolor(curve[i + 1].V, WHITE);
	}
}
void recordCurveFrames(const Curve& curve, VertexRecorder* recorder, float framesize)
{
	Matrix4f T;
	const Vector3f RED(1, 0, 0);
	const Vector3f GREEN(0, 1, 0);
	const Vector3f BLUE(0, 0, 1);
	
	const Vector4f ORGN(0, 0, 0, 1);
	const Vector4f AXISX(framesize, 0, 0, 1);
	const Vector4f AXISY(0, framesize, 0, 1);
	const Vector4f AXISZ(0, 0, framesize, 1);

	for (int i = 0; i < (int)curve.size(); ++i)
	{
		T.setCol(0, Vector4f(curve[i].N, 0));
		T.setCol(1, Vector4f(curve[i].B, 0));
		T.setCol(2, Vector4f(curve[i].T, 0));
		T.setCol(3, Vector4f(curve[i].V, 1));
 
		// Transform orthogonal frames into model space
		Vector4f MORGN  = T * ORGN;
		Vector4f MAXISX = T * AXISX;
		Vector4f MAXISY = T * AXISY;
		Vector4f MAXISZ = T * AXISZ;

		// Record in model space
		recorder->record_poscolor(MORGN.xyz(), RED);
		recorder->record_poscolor(MAXISX.xyz(), RED);

		recorder->record_poscolor(MORGN.xyz(), GREEN);
		recorder->record_poscolor(MAXISY.xyz(), GREEN);

		recorder->record_poscolor(MORGN.xyz(), BLUE);
		recorder->record_poscolor(MAXISZ.xyz(), BLUE);
	}
}

