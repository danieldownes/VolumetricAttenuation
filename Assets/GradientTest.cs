using UnityEngine;
using System.Collections;

public class GradientTest : MonoBehaviour
{

	public GradientColour gradient;

	//public GameObject LinePool;

	public Vector3 startPos;
	public Vector3 gradLegnth;


	void Start()
	{
		gradient = new GradientColour();
		gradient.init();
		gradient.createDefaultHeatMapGradient();	
	}

	void Update()
	{
		float r = 0f;
		float g = 0f;
		float b = 0f;

		for( float n = 1; n >= 0; n = n - 0.01f)
		{
			gradient.getcolourAtValue(n, ref r, ref g, ref b);

			Debug.DrawLine(startPos, startPos + (gradLegnth * n), new Color(r, g, b));
		}
	}

}
