using UnityEngine;
using System;
using System.Collections;
using System.Collections.Generic;

public class HobbyPath : MonoBehaviour
{

	// This class implements a path, which is a list of HobbyPoints
	
	public List<HobbyPoint> points;
	public bool cyclic = true; 	//  Is the path cyclic?
	float curl_begin = 1; 	//  If not, curl parameter at endpoints
	float curl_end = 1;


	void Start()
	{
		Update();
	}

	void Update()
	{
		mp_parse();
		
		solve_angles();
		find_controls();
		
		// Apply tension?
		//foreach(Vector2 pt in points)
		//{
		//p.Append(new Point(pt, 1.0 / tension_, 1.0 / tension_))
		//}
	}

	
	int rangeStart()
	{
		if( cyclic)
			return 0;
		else
			return 1;
	}

	int rangeEnd()
	{
		if( cyclic)
			return points.Count;
		else
			return points.Count - 1;
	}


	
	HobbyPoint getP(int i)
	{
		// Gets the point [i] of the list, but assuming the list is circular and
		//  thus allowing for indexes greater than the list length

		if( i <= -1)
			i = points.Count - 1;

		if( i >= points.Count)
			i = 0;

		return points[i];
	}

	int getI(int i)
	{
		if( i <= -1)
			i = points.Count - 1;
		
		if( i >= points.Count)
			i = 0;

		return i;
	}

	
	// Now some functions from John Hobby and METAFONT book
	// (The good bits)
	
	float f(float theta, float phi)
	{
		// "Velocity" function
		float n = 2 + Mathf.Sqrt(2) * (Mathf.Sin(theta) - Mathf.Sin(phi) / 16) * 
			(Mathf.Sin(phi) - Mathf.Sin(theta) / 16) * 
				(Mathf.Cos(theta) - Mathf.Cos(phi));
		
		float m = 3 * (1 + 0.5f * (Mathf.Sqrt(5) - 1) * Mathf.Cos(theta) + 0.5f * (3 - Mathf.Sqrt(5)) * Mathf.Cos(phi));

		return n / m;
	}
	
	
	Vector2[] control_points(Vector2 z0, Vector2 z1, float theta=0, float phi=0, float alpha=1, float beta=1)
	{
		// Given two points in a path, and the angles of departure and arrival
		// at each one, this function finds the appropriate control points of the
		// Bezier's curve, using John Hobby's algorithm

		Vector2[] ret = new Vector2[2];
		ret[0] = z0 + (Vector2.right * theta).magnitude * (z1 - z0) * f(theta, phi) * alpha;
		ret[1] = z1 - (-Vector2.right * phi ).magnitude * (z1 - z0) * f(phi, theta) * beta;

		return( ret );
	}
	
	void pre_compute_distances_and_angles()
	{
		// This function traverses the path and computes the distance between
		// adjacent points, and the turning angles of the polyline which joins them
		
		for( int i = rangeStart(); i < rangeEnd(); i++ )
		{
			Vector2 v_post = getP(i + 1).z - getP(i).z;
			Vector2 v_ant = getP(i).z - getP(i - 1).z;
			
			// Store the computed values in the Points of the Path
			points[getI(i)].d_ant = v_ant.magnitude;
			points[getI(i)].d_post = v_post.magnitude;
			points[getI(i)].xi = HobbyCurve.arg( HobbyCurve.vec2_quot(v_post, v_ant) );
		}
		
		if( !cyclic)
		{
			// First and last xi are zero
			points[0].xi = getP(-1).xi = 0;
			
			// Also distance to previous and next points are zero for endpoints
			points[0].d_ant = 0;
			points[getI(-1)].d_post = 0;
		}
	}

	
	List<float[]> build_coefficients()
	{
		// This function creates five vectors which are coefficients of a
		// linear system which allows finding the right values of "theta" at
		// each point of the path (being "theta" the angle of departure of the
		// path at each point). The theory is from METAFONT book.
		int c = points.Count;
		float[] A = new float[c];
		float[] B = new float[c]; 
		float[] C = new float[c]; 
		float[] D = new float[c];
		float[] R = new float[c];
		
		pre_compute_distances_and_angles();
		
		if( cyclic)
		{
			//  In this case, first equation doesn't follow the general rule
			float curl = curl_begin;
			float alpha_0 = points[0].alpha;
			float beta_1 = points[1].beta;
			float xi_0 = Mathf.Pow(alpha_0, 2) * curl / Mathf.Pow(beta_1, 2);
			float xi_1 = points[1].xi;
			A[0] = 0;
			B[0] = 0;
			C[0] = (xi_0 * alpha_0) + 3 - beta_1;
			D[0] = (3 - alpha_0) * xi_0 + beta_1;
			R[0] = -D[0] * xi_1;
		}
		
		// Equations 1 to n-1 (or 0 to n for cyclic paths)
		int k;
		for( k = rangeStart(); k < rangeEnd(); k++)
		{
			A[k] = ( getP(k-1).alpha / (Mathf.Pow(getP(k).beta, 2) * getP(k).d_ant));
			B[k] = (3-getP(k-1).alpha) / (Mathf.Pow(getP(k).beta, 2) * getP(k).d_ant);
			        C[k] = ((3-getP(k+1).beta) / (Mathf.Pow(getP(k).alpha, 2) * getP(k).d_post));
			        D[k] = ( getP(k+1).beta / (Mathf.Pow(getP(k).alpha, 2) * getP(k).d_post));
			        R[k] = (-B[k] * getP(k).xi - D[k] * getP(k + 1).xi);
		}
		
		if( !cyclic)
		{
			// The last equation doesn't follow the general form
			int n = R.Length; //  index to generate
			C[k] = 0;
			D[k] = 0;
			float curl = curl_end;
			float beta_n = getP(n).beta;
			float alpha_n_1 = getP(n - 1).alpha;
			float xi_n = Mathf.Pow(beta_n, 2) * curl / Mathf.Pow(alpha_n_1, 2);
			A[k] = ((3-beta_n) * xi_n + alpha_n_1);
			B[k] = (beta_n*xi_n + 3 - alpha_n_1);
			R[k] = (0);
		}
		
		// Prepare Return
		List<float[]> listOut = new List<float[]>();
		listOut.Add(A);
		listOut.Add(B);
		listOut.Add(C);
		listOut.Add(D);
		listOut.Add(R);
		return(listOut);
	}
	
	
	///import numpy as np //  Required to solve the linear equation system
	float[] solve_for_thetas(float[] A, float[] B, float[] C, float[] D, float[] R)
	{
		// This function receives the five vectors created by
		// build_coefficients() and uses them to build a linear system with N
		// unknown (being N the number of points in the path). Solving the system
		// finds the value for theta (departure angle) at each point
		
		int L = R.Length;
		float[,] a = new float[L, L];
		float[] b;

		for(int k = rangeStart(); k < rangeEnd(); k++)
		{
			int prev = k > 0 ? (k - 1) : L - 1;
			int post = (k+1) % L;
			a[k, prev] = A[k];
			a[k, k] = B[k] + C[k];
			a[k, post] = D[k];
			//int[] b = Array.Copy( np.array(R); //TODO: Hard-Copy the array?
		}

		return(ComputeCoefficents(a, R) );  //np.linalg.solve(a, b)
	}
	

	//https://social.msdn.microsoft.com/Forums/en-US/70408584-668d-49a0-b179-fabf101e71e9/solution-of-linear-equations-systems?forum=Vsexpressvcs
	public float[] ComputeCoefficents(float[,] X, float[] Y)
	{
		//float[] Y = new float[X.Length];
	
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

		return Y;
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
		List<float[]> c = build_coefficients();
		float[] x = solve_for_thetas(c[0], c[1], c[2], c[3], c[4]);

		int k;
		for(k = rangeStart(); k < rangeEnd(); k++)
			points[k].theta = x[k];

		for(k = rangeStart(); k < rangeEnd(); k++)
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

		for( int k = rangeStart(); k < rangeEnd(); k++)
		{
			Vector2 z0 = getP(k).z;
			Vector2 z1 = getP(k + 1).z;
			float theta = getP(k).theta;
			float phi = getP(k + 1).phi;
			float alpha = getP(k).alpha;
			float beta = getP(k + 1).beta;
			Vector2[] uv = control_points(z0, z1, theta, phi, alpha, beta);
			points[getI(k)].u_right = uv[0];
			points[getI(k+1)].v_left = uv[1];
		}
	}

		

	void mp_parse()
	{
		// This function receives a string which contains a path in metapost syntax,
		// and returns a Path object which stores the same path in the structure
		// required to compute the control points.
		// The path should only contain explicit coordinates and numbers.
		// Currently only "curl" and "tension" keywords are understood. Direction
		// options are ignored.
		float default_tension = 1;
		float default_curl = 1;

		float alpha;
		float beta;

		curl_begin = default_curl;
		curl_end = default_curl;
		alpha = beta = (1.0f / default_tension);

		int k = 0;
		foreach(HobbyPoint p in points)
		{
			if( p.tension != 0)
			{
				alpha = 1.0f / p.tension;
			}
			if( points.Count > 3) ///len(aux)>3:
				beta = 1.0f / p.tension; ///float(aux[3])
			else
				beta = alpha;

			if( p.curl != 0)
			{
				if(k == 0)
					curl_begin = p.curl; ///float(aux[1])
				else
					curl_end = p.curl; ///float(aux[1])
			}

			//  Update tensions
			points[getI(k-1)].alpha = alpha;
			points[getI(k)].beta = beta;
		}

	}

}
