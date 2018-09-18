#include "surf.h"
#include "vertexrecorder.h"
using namespace std;

const float c_pi = 3.14159265358979323846f;

namespace
{
    
    // We're only implenting swept surfaces where the profile curve is
    // flat on the xy-plane.  This is a check function.
    static bool checkFlat(const Curve &profile)
    {
        for (unsigned i=0; i<profile.size(); i++)
            if (profile[i].V[2] != 0.0 ||
                profile[i].T[2] != 0.0 ||
                profile[i].N[2] != 0.0)
                return false;
    
        return true;
    }
}

// DEBUG HELPER
Surface quad() { 
	Surface ret;
	ret.VV.push_back(Vector3f(-1, -1, 0));
	ret.VV.push_back(Vector3f(+1, -1, 0));
	ret.VV.push_back(Vector3f(+1, +1, 0));
	ret.VV.push_back(Vector3f(-1, +1, 0));

	ret.VN.push_back(Vector3f(0, 0, 1));
	ret.VN.push_back(Vector3f(0, 0, 1));
	ret.VN.push_back(Vector3f(0, 0, 1));
	ret.VN.push_back(Vector3f(0, 0, 1));

	ret.VF.push_back(Tup3u(0, 1, 2));
	ret.VF.push_back(Tup3u(0, 2, 3));
	return ret;
}

Surface makeSurfRev(const Curve &profile, unsigned steps)
{
    Surface surface;
//	surface = quad();
    
    if (!checkFlat(profile))
    {
        cerr << "surfRev profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    // TODO: Here you should build the surface.  See surf.h for details.
    //steps = 10;
    float step_size = 2*c_pi/steps;
    
    for(int i = 0; i < steps; i++)
    {
        float theta = i*step_size;
//        std::cout << "theta " << theta << std::endl;
        Matrix3f R( cos(theta), 0, sin(theta),
                    0         , 1, 0,
                   -sin(theta), 0, cos(theta));

        Matrix3f NR = R.inverse();
        NR.transpose();
        
        for (int j = 0; j < (int)profile.size(); j++)
        {
            Vector3f vertex = profile[j].V;
            vertex = R*vertex;
            surface.VV.push_back(vertex);
            
            Vector3f normal = profile[j].N;
            normal = -1*NR*normal;
            surface.VN.push_back(normal);
            
//            std::cout << vertex[0] << " " << vertex[1] << " " << vertex[2] << std::endl;
        }
        
//        std::cout << "vertices " << surface.VV.size() << std::endl;
        
    }
    
    
    int n_vertices = profile.size();
    //steps = 3;
    //n_vertices = 6;
    for(int k = 0; k < steps - 1; k++)
    {
        for(int l = 1; l < n_vertices; l++)
        {
            int m = n_vertices*k;
            surface.VF.push_back(Tup3u(m + l, m + l + n_vertices, m + l + n_vertices - 1));
//            std::cout << m + l << " " << m + l + n_vertices - 1 << " " << m + l - 1 << std::endl;
//            std::cout << m + l << " " << m + l + n_vertices << " " << m + l + n_vertices - 1 << std::endl;
            surface.VF.push_back(Tup3u(m + l, m + l + n_vertices - 1, m + l - 1));
        }
    }
    //the last case: k = steps-1
    for(int l = 1; l < n_vertices; l++)
    {
        int m = n_vertices*(steps-1);
        surface.VF.push_back(Tup3u(m+l, l , l - 1));
//        std::cout << m+l << " " << l << " " << l - 1 << std::endl;
        surface.VF.push_back(Tup3u(m+l, l - 1, m + l - 1));
//        std::cout << m+l << " " << l-1 << " " << m + l - 1 << std::endl;
//        std::cout << "vertices " << surface.VV.size() << std::endl;
    }
    
    std::cout << "steps * n_vertices " << steps*n_vertices << " total vertices " << surface.VV.size() << std::endl;
//    cerr << "\t>>> makeSurfRev called (but not implemented).\n\t>>> Returning empty surface." << endl;
 
    return surface;
//    return Surface();
}

void addFaces(const Curve &profile, int steps, Surface* surface)
{
    int n_vertices = profile.size();
    //steps = 3;
    //n_vertices = 6;
    for(int k = 0; k < steps - 1; k++)
    {
        for(int l = 1; l < n_vertices; l++)
        {
            int m = n_vertices*k;
            surface->VF.push_back(Tup3u(m + l, m + l + n_vertices, m + l + n_vertices - 1));
            //            std::cout << m + l << " " << m + l + n_vertices - 1 << " " << m + l - 1 << std::endl;
            //            std::cout << m + l << " " << m + l + n_vertices << " " << m + l + n_vertices - 1 << std::endl;
            surface->VF.push_back(Tup3u(m + l, m + l + n_vertices - 1, m + l - 1));
        }
    }
    //the last case: k = steps-1
    for(int l = 1; l < n_vertices; l++)
    {
        int m = n_vertices*(steps-1);
        surface->VF.push_back(Tup3u(m+l, l , l - 1));
        //        std::cout << m+l << " " << l << " " << l - 1 << std::endl;
        surface->VF.push_back(Tup3u(m+l, l - 1, m + l - 1));
        //        std::cout << m+l << " " << l-1 << " " << m + l - 1 << std::endl;
        //        std::cout << "vertices " << surface.VV.size() << std::endl;
    }
}

Surface makeGenCyl(const Curve &profile, const Curve &sweep )
{
    Surface surface;
//	surface = quad();

    if (!checkFlat(profile))
    {
        cerr << "genCyl profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    int steps = sweep.size();
    for(int i = 0; i < steps; i++)
    {

        Vector4f sweepN(sweep[i].N, 0);
        Vector4f sweepB(sweep[i].B, 0);
        Vector4f sweepT(sweep[i].T, 0);
        Vector4f sweepV(sweep[i].V, 1);
        Matrix4f M(sweepN, sweepB, sweepT, sweepV);
        
        Matrix3f NR = M.getSubmatrix3x3(0, 0).inverse();
        NR.transpose();
        
        for (int j = 0; j < (int)profile.size(); j++)
        {
            Vector4f profile_vertex = Vector4f(profile[j].V, 1);
            Vector3f vertex = (M*profile_vertex).xyz();
            surface.VV.push_back(vertex);
            
            Vector3f normal = profile[j].N;
            normal = -1*NR*normal;
            surface.VN.push_back(normal);
        }
    }
    
    addFaces(profile, steps, &surface);
    
    cerr << "\t>>> makeGenCyl called (but not implemented).\n\t>>> Returning empty surface." <<endl;
    
    return surface;
}

void recordSurface(const Surface &surface, VertexRecorder* recorder) {
	const Vector3f WIRECOLOR(0.4f, 0.4f, 0.4f);
    for (int i=0; i<(int)surface.VF.size(); i++)
    {
		recorder->record(surface.VV[surface.VF[i][0]], surface.VN[surface.VF[i][0]], WIRECOLOR);
		recorder->record(surface.VV[surface.VF[i][1]], surface.VN[surface.VF[i][1]], WIRECOLOR);
		recorder->record(surface.VV[surface.VF[i][2]], surface.VN[surface.VF[i][2]], WIRECOLOR);
    }
}

void recordNormals(const Surface &surface, VertexRecorder* recorder, float len)
{
	const Vector3f NORMALCOLOR(0, 1, 1);
    for (int i=0; i<(int)surface.VV.size(); i++)
    {
		recorder->record_poscolor(surface.VV[i], NORMALCOLOR);
		recorder->record_poscolor(surface.VV[i] + surface.VN[i] * len, NORMALCOLOR);
    }
}

void outputObjFile(ostream &out, const Surface &surface)
{
    
    for (int i=0; i<(int)surface.VV.size(); i++)
        out << "v  "
            << surface.VV[i][0] << " "
            << surface.VV[i][1] << " "
            << surface.VV[i][2] << endl;

    for (int i=0; i<(int)surface.VN.size(); i++)
        out << "vn "
            << surface.VN[i][0] << " "
            << surface.VN[i][1] << " "
            << surface.VN[i][2] << endl;

    out << "vt  0 0 0" << endl;
    
    for (int i=0; i<(int)surface.VF.size(); i++)
    {
        out << "f  ";
        for (unsigned j=0; j<3; j++)
        {
            unsigned a = surface.VF[i][j]+1;
            out << a << "/" << "1" << "/" << a << " ";
        }
        out << endl;
    }
}
