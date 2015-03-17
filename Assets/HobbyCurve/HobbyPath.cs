using UnityEngine;
using System;
using System.Collections;
using System.Collections.Generic;

public class HobbyPath : MonoBehaviour
{

	// This class implements a path, which is a list of HobbyPoints
	
	public List<HobbyPoint> points;
	bool cyclic = true; 	//  Is the path cyclic?
	int curl_begin = 1; 	//  If not, curl parameter at endpoints
	int curl_end = 1;
	
	public void init()
	{
		init(1, true, 1, 1);
	}
	
	public void init(float tension_, bool cyclic_, int curl_begin, int curl_end_)
	{
		points = List<HobbyPoint>();
		
		cyclic = cyclic_;
		curl_begin = curl_begin;
		curl_end = curl_end_;
		
		
		// Apply tension?
		//foreach(Vector2 pt in points)
		//{
		//p.Append(new Point(pt, 1.0 / tension_, 1.0 / tension_))
		//}
	}
	
	
	int[] range()
	{
		// Returns the range of the indexes of the points to be solved.
		// This range is the whole length of p for cyclic paths, but excludes
		// the first and last points for non-cyclic paths
		
		if( cyclic)
			return range(points.Count);
		else
			return range(1, points.Count - 1);
	}
	
	// The following functions allow to use a Path object like an array
	// so that, if x = Path(...), you can do len(x) and x[i]
	void append(List<Point> data)
	{
		points.AddRange(data);
	}
	
	
	void len()
	{
		return points.Count;
	}
	
	
	void getitem(int i)
	{
		// Gets the point [i] of the list, but assuming the list is circular and
		//  thus allowing for indexes greater than the list length
		i %= points.Count;
		return points[i];
	}

	
	// Now some functions from John Hobby and METAFONT book
	// (The good bits)
	
	void f(float theta, float phi)
	{
		// "Velocity" function
		float n = 2 + Mathf.Sqrt(2) * (Mathf.Sin(theta) - Mathf.Sin(phi) / 16) * 
			(Mathf.Sin(phi) - Mathf.Sin(theta) / 16) * 
				(Mathf.Cos(theta) - Mathf.Cos(phi));
		
		float m = 3 * (1 + 0.5f * (Mathf.Sqrt(5) - 1) * Mathf.Cos(theta) + 0.5f * (3 - Mathf.Sqrt(5)) * Mathf.Cos(phi));
		return n / m;
	}
	
	
	Vector2 control_points(float z0, float z1, float theta=0, float phi=0, float alpha=1, float  beta=1)
	{
		// Given two points in a path, and the angles of departure and arrival
		// at each one, this function finds the appropriate control points of the
		// Bezier's curve, using John Hobby's algorithm
		
		Vector2 i = Vector2(0, 1);
		float u = z0 + Mathf.Exp(i * theta) * (z1 - z0) * f(theta, phi) * alpha;
		float v = z1 - Mathf.Exp(-i * phi) * (z1 - z0) * f(phi, theta) * beta;
		
		return( new Vector2(u, v) );
	}
	
	
	void pre_compute_distances_and_angles()
	{
		// This function traverses the path and computes the distance between
		// adjacent points, and the turning angles of the polyline which joins them
		
		for( int i = 0; i < points.Count; i++ )
		{
			float v_post = points[i + 1].z - points[i].z;
			float v_ant = points[i].z - points[i - 1].z;
			
			// Store the computed values in the Points of the Path
			points[i].d_ant = Mathf.Abs(v_ant);
			points[i].d_post = Mathf.Abs(v_post);
			points[i].xi = HobbyCurve.arg(v_post / v_ant);
		}
		
		if( !cyclic)
		{
			// First and last xi are zero
			points[0].xi = points[-1].xi = 0;
			
			// Also distance to previous and next points are zero for endpoints
			points[0].d_ant = 0;
			points[-1].d_post = 0;
		}
	}
	
	List<float[]> build_coefficients()
	{
		// This function creates five vectors which are coefficients of a
		// linear system which allows finding the right values of "theta" at
		// each point of the path (being "theta" the angle of departure of the
		// path at each point). The theory is from METAFONT book.
		
		float[] A;
		float[] B; 
		float[] C; 
		float[] D;
		float[] R;
		
		pre_compute_distances_and_angles();
		
		if( cyclic)
		{
			//  In this case, first equation doesn't follow the general rule
			int curl = curl_begin;
			float alpha_0 = points[0].alpha;
			float beta_1 = points[1].beta;
			float xi_0 = (alpha_0 ** 2) * curl / (beta_1 ** 2);
			float xi_1 = points[1].xi;
			A[0] = 0;
			B[0] = 0;
			C[0] = (xi_0 * alpha_0) + 3 - beta_1;
			D[0] = (3 - alpha_0) * xi_0 + beta_1;
			R[0] = -D[0] * xi_1;
		}
		
		// Equations 1 to n-1 (or 0 to n for cyclic paths)
		int k;
		for( k = 1; k < points.Count; k++) //TODO: Replace with range()
		{
			A[k] = ( points[k-1].alpha / ((points[k].beta**2) * points[k].d_ant));
			B[k] = (3-points[k-1].alpha) / ((points[k].beta**2) * points[k].d_ant);
			C[k] = ((3-points[k+1].beta) / ((points[k].alpha**2) * points[k].d_post));
			D[k] = ( points[k+1].beta / ((points[k].alpha**2) * points[k].d_post));
			R[k] = (-B[k] * points[k].xi - D[k] * points[k+1].xi);
		}
		
		if( !cyclic)
		{
			// The last equation doesn't follow the general form
			int n = R.Length; //  index to generate
			C[k] = 0;
			D[k] = 0;
			int curl = curl_end;
			float beta_n = points[n].beta;
			float alpha_n_1 = points[n - 1].alpha;
			float xi_n = (beta_n**2) * curl / (alpha_n_1**2); 	// TODO: ** ?????
			A[k]((3-beta_n) * xi_n + alpha_n_1);
			B[k](beta_n*xi_n + 3 - alpha_n_1);
			R[k](0);
		}
		
		// Prepare Return
		List<float[]> listOut = new List<float[]>();
		listOut.Add(A, B, C, D, R);
		return(listOut);
	}
	
	
	///import numpy as np //  Required to solve the linear equation system
	void solve_for_thetas(float[] A, float[] B, float[] C, float[] D, float[] R)
	{
		// This function receives the five vectors created by
		// build_coefficients() and uses them to build a linear system with N
		// unknown (being N the number of points in the path). Solving the system
		// finds the value for theta (departure angle) at each point
		
		int L = R.Length;
		int[,] a = new int[L, L];

		for(int k = 0; k < L; k++)
		{
			int prev = (k - 1) % L;
			int post = (k + 1) % L;
			a[k][prev] = A[k];
			a[k][k] = B[k] + C[k];
			a[k][post] = D[k];
			int[] b = np.array(R); //TODO: Hard-Copy the array?
		}	

		return(ComputeCoefficents(a, b) );  //np.linalg.solve(a, b)
	}
	

	//https://social.msdn.microsoft.com/Forums/en-US/70408584-668d-49a0-b179-fabf101e71e9/solution-of-linear-equations-systems?forum=Vsexpressvcs
	public void ComputeCoefficents(float[,] X, float[] Y)
	{
		  int I, J, K, K1, N;
		  N = Y.Length;
		  for (K = 0; K < N; K++)
		  {
			K1 = K + 1;
			for (I = K; I < N; I++)
			{
			  if (X[I, K] != 0)
			  {
				for (J = K1; J < N; J++)
				{
				  X[I, J] /= X[I, K];
				}
				Y[I] /= X[I, K];
			  }
			}
			for (I = K1; I < N; I++)
			{
			  if (X[I, K] != 0)
			  {
				for (J = K1; J < N; J++)
				{
				  X[I, J] -= X[K, J];
				}
				Y[I] -= Y[K];
			  }
			}
		  }
		  for (I = N - 2; I >= 0; I--)
		  {
			for (J = N - 1; J >= I + 1; J--)
			{
			  Y[I] -= X[I, J] * Y[J];
			}
		  }
		}
	}
	
	
	
	void solve_angles()
	{
		// This function receives a path in which each point is "open", i.e. it
		// does not specify any direction of departure or arrival at each node,
		// and finds these directions in such a way which minimizes "mock
		// curvature". The theory is from METAFONT book.
		
		// Basically it solves
		// a linear system which finds all departure angles (theta), and from
		// these and the turning angles at each point, the arrival angles (phi)
		// can be obtained, since theta + phi + xi = 0 at each knot
		float x = solve_for_thetas(*build_coefficients());  //TODO: is * needed?
		
		for(int k = 0; k < points.Count; k++)
			points[k].theta = x[k];

		foreach(int k in range(L))
			points[k].phi = - points[k].theta - points[k].xi;
	}
	
	
	void find_controls()
	{
		// This function receives a path in which, for each point, the values
		// of theta and phi (leave and enter directions) are known, either because
		// they were previously stored in the structure, or because it was
		// computed by function solve_angles(). From this path description
		// this function computes the control points for each knot and stores
		// it in the path. After this, it is possible to print path to get
		// a string suitable to be feed to tikz.
		
		///r = []
		for( int k; k < points.Count; k++)  //TODO: range() ??
		{
			float z0 = points[k].z;
			float z1 = points[k + 1].z;
			float theta = points[k].theta;
			float phi = points[k + 1].phi;
			float alpha = points[k].alpha;
			float beta = points[k + 1].beta;
			Vector2 uv = control_points(z0, z1, theta, phi, alpha, beta);
			points[k].u_right = uv.x;
			points[k + 1].v_left = uv.y;
		}
	}


	
	/*
	void str__(self)
	{
		// String serialization
		// The printable representation of the object is one suitable for
		// feeding it into tikz, producing the same figure than in metapost
		
		///r = []
		///L = len(self.p)
		last = 1;
		if(self.cyclic)
			last = 0;
		
		for( k in range(L-last) )
		{
			post = (k + 1) % L;
			z = self.p[k].z'
			u = self.p[k].u_right;
			v = self.p[post].v_left;
			r.append("(%.4f, %.4f) .. controls (%.5f, %.5f) and (%.5f, %.5f)" %\
				(z.real, z.imag, u.real, u.imag, v.real, v.imag))
			if(self.cyclic:
				last_z = self.p[0].z
		else
			last_z = self.p[-1].z
		r.append("(%.4f, %.4f)" % (last_z.real, last_z.imag))
		
		return "..".join(r);
	}
	*/
	
	
	/*
	void repr(self)
	{
		// Dumps internal parameters, for debugging purposes
		///r = ["Path information"]
		///r.append("Cyclic=%s, curl_begin=%s, curl_end=%s" % (self.cyclic,
		///self.curl_begin, self.curl_end))
		///for pt in self.p:
		///r.append(str(pt))
		///return "\n".join(r)
	}
	*/

	
	//void mp_to_tikz(path, command=None, options=None)
	//{
	/*
		Utility receives a string containing a metapost path
		and uses all the above to generate the tikz version with explicit
		control points.
		It does not make a full parsing of the metapost path. Currently it is
		not possible to specify directions nor tensions at knots. It uses
		default tension = 1, default curl =1 for both ends in non-cyclic paths
		and computes the optimal angles at each knot. It does admit however
		cyclic and non-cyclic paths.
		To summarize, the only allowed syntax is z0 .. z1 .. z2, where z0, z1,
		etc are explicit coordinates such as (0,0) .. (1,0) etc.. And
		optionally the path can ends with the literal "cycle".
		*/
	///tension = 1;
	///curl = 1;
	///if options:
	///opt = []
	///for o in options.split(","):
	///o=o.strip()
	///if o.startswith("tension"):
	///tension = float(o.split("=")[1])
	///elif o.startswith("curl"):
	///curl = float(o.split("=")[1])
	///else:
	///opt.append(o)
	///options = ",".join(opt)
	///new_path = mp_parse(path, default_tension = tension, default_curl = curl)
	/////  print repr(new_path)
	///solve_angles(new_path)
	///find_controls(new_path)
	///if command==None:
	///command="draw"
	///if options==None:
	///options = ""
	///else:
	///options = "[%s]" % options
	
	//return "\\%s%s %s;" % (command, options, str(new_path))
	//}
	
	/* // Needs converting?
	void mp_parse(mppath, default_tension = 1, default_curl = 1)
	{
		// This function receives a string which contains a path in metapost syntax,
		// and returns a Path object which stores the same path in the structure
		// required to compute the control points.
		// The path should only contain explicit coordinates and numbers.
		// Currently only "curl" and "tension" keywords are understood. Direction
		// options are ignored.

		if mppath.endswith(";"): //  Remove last semicolon
		mppath=mppath[:-1]
		pts = mppath.split("..") //  obtain points
		pts = [p.strip() for p in pts] //  remove extra spaces

		if pts[-1] == "cycle":
			is_cyclic = True
			pts=pts[:-1] //  Remove this last keyword
		else:
			is_cyclic = False
			path = Path([], cyclic=is_cyclic)
			path.curl_begin = default_curl
			path.curl_end = default_curl
			alpha = beta = 1.0/default_tension
			k=0
			for p in pts:
			if p.startswith("tension"):
				aux = p.split()
				alpha = 1.0/float(aux[1])	//TODO: nesting ????
			if len(aux)>3:
			beta = 1.0/float(aux[3])
			else:
			beta = alpha
			else:
			aux = p.split("{") //  Extra options at the point
			p = aux[0].strip()
			if p.startswith("curl"):
			if k==0:
			path.curl_begin=float(aux[1])
			else:
			path.curl_end = float(aux[1])
			elif p.startswith("dir"):
			//  Ignored by now
			pass
		
		path.append(Point(eval(p))) //  store the pair of coordinates
		//  Update tensions
		path[k-1].alpha = alpha
		path[k].beta = beta
	}
	*/
}
