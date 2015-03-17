using UnityEngine;
using System;
using System.Collections;
using System.Collections.Generic;

public class HobbyPoint : MonoBehaviour
{

	// This class implements the coordinates of a knot, and all kind of auxiliar parameters to compute a smooth path passing through it
	public Vector2 z = Vector2(0, 0); 		//  Point coordinates
	public float alpha = 1; 				//  Tension at point (1 by default)
	public float beta = 1;
	public float theta = 0; 				//  Angle at which the path leaves
	public float phi = 0; 					//  Angle at which the path enters
	public float xi = 0; 					//  angle turned by the polyline at this point
	public Vector2 v_left = Vector2(0,0); 	//  Control points of the Bezier curve at this point
	public Vector2 u_right = Vector2(0,0); 	//  (to be computed later)
	public float d_ant = 0; 				//  Distance to previous point in the path
	public float d_post = 0; 				//  Distance to next point in the path
	
	
	// May require alternative function for optional parameters
	// Coordinates can be given as a complex number
	public void init(Vector2 z_)
	{
		init(z_, 1, 1, Vector2(0, 0),Vector2(0, 0) );
	}
	public void init(Vector2 z_, float alpha_, float beta_, Vector2 v_, Vector2 u_)
	{
		z = z_;
		alpha = alpha_;
		beta = beta_;
		v_left = v_;
		u_right = u_;
		d_ant = 0;
		d_post = 0;
		xi = 0;
	}
	
	
	/*
	void  __str__(self)
	{
		// Creates a printable representation of this object, for
		// debugging purposes return
		///z = (%.3f, %.3f) alpha = %.2f beta = %.2f theta=%.2f phi=%.2f
		///[v=(%.2f, %.2f) u=(%.2f, %.2f) d_ant=%.2f d_post=%.2f xi=%.2f]""" % \
		///(self.z.real, self.z.imag, self.alpha, self.beta,
		///degrees(self.theta), degrees(self.phi),
		///self.v_left.real, self.v_left.imag, self.u_right.real,
		///self.u_right.imag, self.d_ant, self.d_post, degrees(self.xi))
	}
	*/
}
