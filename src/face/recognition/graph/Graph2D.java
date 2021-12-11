package face.recognition.graph;

import java.awt.Graphics;

import javax.swing.JFrame;

public class Graph2D extends JFrame
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	int [] xPoints = null;
	int [] yPoints = null;
	
	public Graph2D(String title)
	{
		this.setTitle(title);
	}
	public void setPoints(int [] xPoints, int [] yPoints)
	{
		this.xPoints = xPoints;
		this.yPoints = yPoints;
	}
	
	public void paint(Graphics gDrawer)
	{
		int panelW = this.getWidth() - 10;
		int panelH = this.getHeight() - 10;
		
		int maxX = 0;
		int maxY = 0;
		
		for(int i = 0; i < xPoints.length; i++)
		{
			if(maxX < xPoints[i])
			{
				maxX = xPoints[i];
			}
		}
		
		for(int i = 0; i < yPoints.length; i++)
		{
			if(maxY < yPoints[i])
			{
				maxY = yPoints[i];
			}
		}
		
		double xUnit = (double)(panelW - 20) / maxX;
		double yUnit = (double)(panelH - 20) / maxY;
		
		int startGraphX = 5;
		int startGraphY = panelH - 5;

		try
		{
			for(int i = 1; i < xPoints.length; i++)
			{
				int startLineX = startGraphX + (int)(xPoints[i - 1] * xUnit);
				int startLineY = startGraphY - (int)(yPoints[i - 1] * yUnit);
				
				int endLineX = startGraphX + (int)(xPoints[i] * xUnit);
				int endLineY = startGraphY - (int)(yPoints[i] * yUnit);
				
				gDrawer.drawLine(startLineX, startLineY, endLineX, endLineY);
			}
		}
		catch(Exception ex) 
		{
			System.out.println(ex.getMessage());
		}
	}
}
