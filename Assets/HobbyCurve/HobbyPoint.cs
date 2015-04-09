using UnityEngine;
using System;
using System.Collections;
using System.Collections.Generic;

public class HobbyPoint : MonoBehaviour
{

	// This class implements the coordinates of a knot, and all kind of auxiliar parameters to compute a smooth path passing through it
	public Vector2 z = new Vector2(0, 0); 		//  Point coordinates
	public float alpha = 1; 				//  Tension at point (1 by default)
	public float beta = 1;
	public float theta = 0; 				//  Angle at which the path leaves
	public float phi = 0; 					//  Angle at which the path enters
	public float xi = 0; 					//  angle turned by the polyline at this point
	public Vector2 v_left = new Vector2(0,0); 	//  Control points of the Bezier curve at this point
	public Vector2 u_right = new Vector2(0,0); 	//  (to be computed later)
	public float d_ant = 0; 				//  Distance to previous point in the path
	public float d_post = 0; 				//  Distance to next point in the path


	public float tension = 0.5f;
	public float curl = 0.5f;

	public Transform cp_left;
	public Transform cp_right;


	void Start()
	{
		// Grab control point transforms
		foreach (Transform t in transform)
		{
			if(t.name == "c0")
				cp_left = t;
			else if(t.name == "c1")
				cp_right = t;
		}
		
		Update();
		
		d_ant = 0;
		d_post = 0;
		xi = 0;
	}

	void Update()
	{
		z.x = transform.position.x;
		z.y = transform.position.y;

		//cp_left.position = new Vector3(v_left.x, v_left.y, 0);
		//cp_right.position = new Vector3(u_right.x, u_right.y, 0);
	}


}
