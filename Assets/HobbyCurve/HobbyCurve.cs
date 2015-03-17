using UnityEngine;
using System;
using System.Collections;
using System.Collections.Generic;

// Implementation of John Hobby's algorithm

/*
 Find a smooth curve which passes through a series of points.
 The algorithm is used in METAFONT and MetaPost, but the source code
  of these programs is hard to read. I tried to implement it in a more
  modern way, which makes the algorithm more understandable and perhaps portable
  to other languages

  For the second case, the use is:
 
  $ python mp2tikz.py <metapost path> <options>
 
  Where:
  <metapost path> is a path using metapost syntax with the following restrictions:
  * All points have to be explicit (no variables or expressions)
  * All joins have to be "curved" ( .. operator)
  * Options in curly braces next to the nodes are ignored, except
  for {curl X} at end points
  * tension can be specified using metapost syntax
  * "cycle" as end point denotes a cyclic path, as in metapost
  Examples:
  (0,0) .. (60,40) .. (40,90) .. (10,70) .. (30,50) .. cycle
  (0,0) .. (60,40) .. (40,90) .. (10,70) .. (30,50)
  (0,0){curl 10} .. (60,40) .. (40,90) .. (10,70) .. (30,50)
  (0,0) .. (60,40) .. (40,90) .. tension 3 .. (10,70) .. (30,50) .. cycle
  (0,0) .. (60,40) .. (40,90) .. tension 1 and 3 .. (10,70) .. (30,50) .. cycle
 
  <options> can be:
  tension = X. The given tension is applied to all segments in the path by default
  (but tension given at specific points override this setting at those points)
  curl = X. The given curl is applied by default to both ends of the open path
  (but curl given at specific endings override this setting at that point)
  any other options are considered tikz options.
 
  The script prints in standard output a tikz command which draws the given path
  using the given options. In this path all control points are explicit, as computed
  by the string using Hobby's algorithm.
 
  For example:
 
  $ python mp2tikz.py "(0,0) .. (10,10) .. (20,0) .. (10, -10) .. cycle" "tension =3, blue"
 
  Would produce
  \draw[blue] (0.0000, 0.0000) .. controls (-0.00000, 1.84095) and (8.15905, 10.00000)..
  (10.0000, 10.0000) .. controls (11.84095, 10.00000) and (20.00000, 1.84095)..
  (20.0000, 0.0000) .. controls (20.00000, -1.84095) and (11.84095, -10.00000)..
  (10.0000, -10.0000) .. controls (8.15905, -10.00000) and (0.00000, -1.84095)..(0.0000, 0.0000);



 Conversion notes:
	**  is Mathf.Pow
 
	numpy.linalg.solve(a, b)
	
	complex type == Vector2
 
 
   Unsure: 
   radians() in Phython is Deg2Rad?
 
 
  range(list) returns array of indexes in the list
  a = np.zeros((L, L)); // Double brackets == 2D array
  solve_for_thetas(*build_coefficients(path))    * is for???
  r = []  //  before for loops ??
  points[k].phi = -    is -= ??
  
 */


public class HobbyCurve : MonoBehaviour
{
	// Helper functions..

	static float arg(Vector2 z)
	{
		return Mathf.Atan2(z.x, z.y);
	}

	static Vector2 direc(float angle)
	{
		// Given an angle in degrees, returns a complex with modulo 1 and the given phase
		float phi = Mathf.Deg2Rad(angle); //TODO: radians() in Python is Deg2Rad?
		return Vector2(Mathf.Cos(phi), Mathf.Sin(phi));
	}

	static void direc_rad(float angle)
	{
		// Given an angle in radians, returns a complex with modulo 1 and the given phase
		return Vector2(Mathf.Cos(angle), Mathf.Sin(angle));
	}
}