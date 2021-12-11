package face.recognition.pca;

/*-----------------------------------------------------------------------------
File         : MainFrame
 Author      : Engin DURMAZ
 Date        : 23.04.2019
 Description : Class which manages graphical user interface
              Get inputs from user like training path and test path
              Call agentPCA to train images
              Call agentPCA to find nearest label
-----------------------------------------------------------------------------*/

import java.awt.Color;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.SystemColor;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.JTextPane;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.LineBorder;

import face.recognition.graph.Graph2D;

public class MainFrame extends JFrame
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	
	private JPanel mainContentPane;
	private JPanel pnlInputs = new JPanel();
	private JTextPane txtDebug = new JTextPane();

	private JTextField txtTrainingImagePath = null;// path for training images
	private JTextField txtTestImagePath = null; // path for single test image
	private JTextField txtTestAllPath = null; // path for all test images folder
	private JSpinner spinFeatureDimension= null;	

	private int MAIN_FRAME_XPOS = 50;
	private int MAIN_FRAME_YPOS = 50;
	private int MAIN_FRAME_WIDTH = 800;
	private int MAIN_FRAME_HEIGHT = 600;
	

	private int MARGIN_FRAME = 2;
	private int INPUT_PANEL_HEIGHT = 150;
	private int INPUT_PANEL_WIDTH = MAIN_FRAME_WIDTH - 2 * MARGIN_FRAME;
	
	private int DEBUG_PANEL_YPOS = INPUT_PANEL_HEIGHT + 2 * MARGIN_FRAME;
	private int DEBUG_PANEL_HEIGHT = MAIN_FRAME_HEIGHT - INPUT_PANEL_HEIGHT - 30;
	private int DEBUG_PANEL_WIDTH = MAIN_FRAME_WIDTH - 2 * MARGIN_FRAME;

	Color MAIN_BACKGROUND_COLOR = new Color(100, 100, 150);
	Color DEBUG_BACKGROUND_COLOR = new Color(240, 230, 140);
	AgentPCA trainer = null;
		
	/*
	 *  if it is wanted to be run from 1 to 100 dimension, enable this variable
	 *  isForAllDimensions is enabled/disabled by checkbox on screen
	 */
	private boolean isForAllDimensions = false;
	
	/*
	 * if debug enabled print all processes to textbox on the form
	 * if not enabled, do not update text box; it is disabled when running all dimensions for all
	 * because of faster results 
	 */
	private boolean isDebugEnabled = true;
	
	public MainFrame() 
	{
		setTitle("PCA Recognition");
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(MAIN_FRAME_XPOS, MAIN_FRAME_YPOS, MAIN_FRAME_WIDTH, MAIN_FRAME_HEIGHT);
		mainContentPane = new JPanel();
		mainContentPane.setBackground(SystemColor.inactiveCaption);
		mainContentPane.setBorder(new LineBorder(new Color(0, 0, 0)));
		setContentPane(mainContentPane);
		mainContentPane.setLayout(null);
				

		pnlInputs.setBounds(MARGIN_FRAME, MARGIN_FRAME, INPUT_PANEL_WIDTH, INPUT_PANEL_HEIGHT);
		mainContentPane.setBackground(MAIN_BACKGROUND_COLOR);
		mainContentPane.add(pnlInputs);
		pnlInputs.setLayout(null);
		
		JLabel lblFeatureDimension = new JLabel("Feature Dimension");
		lblFeatureDimension.setBounds(5, 5, 150, 25);
		pnlInputs.add(lblFeatureDimension);
		
		spinFeatureDimension = new JSpinner(new SpinnerNumberModel(3, 1, 100, 1));
		spinFeatureDimension.setBounds(200, 5, 200, 25);
		pnlInputs.add(spinFeatureDimension);
		
		JCheckBox chkForallDimensions = new JCheckBox("For all dimensions");
		chkForallDimensions.setBounds(450, 5, 150, 25);
		pnlInputs.add(chkForallDimensions);
		chkForallDimensions.addActionListener
		(
				new ActionListener() {
					public void actionPerformed(ActionEvent e) 
					{
						isForAllDimensions = chkForallDimensions.isSelected();
					}
				}
		);
		
		JLabel lblTrainingImagePath = new JLabel("Training Image Path");
		lblTrainingImagePath.setBounds(5, 35, 150, 25);
		pnlInputs.add(lblTrainingImagePath);
		
		txtTrainingImagePath = new JTextField(1085);
		txtTrainingImagePath.setBounds(200, 35, 200, 25);
		pnlInputs.add(txtTrainingImagePath);
		
		JButton btnTraining = new JButton("Train");
		btnTraining.setBounds(450, 35, 100, 25);
		pnlInputs.add(btnTraining);
		btnTraining.addActionListener
		(
				new ActionListener() {
					public void actionPerformed(ActionEvent e) 
					{
						onTrainingClicked();
					}
				}
		);
		
		JLabel lblTestImagePath = new JLabel("Test Image Path");
		lblTestImagePath.setBounds(5, 65, 150, 25);
		pnlInputs.add(lblTestImagePath);
		
		txtTestImagePath = new JTextField(1085);
		txtTestImagePath.setBounds(200, 65, 200, 25);
		pnlInputs.add(txtTestImagePath);
		
		JButton btnTest = new JButton("Test");
		btnTest.setBounds(450, 65, 100, 25);
		pnlInputs.add(btnTest);
		btnTest.addActionListener
		(
				new ActionListener() {
					public void actionPerformed(ActionEvent e) 
					{
						onTestClicked();
					}
				}
		);

		txtTestAllPath = new JTextField(1085);
		txtTestAllPath.setBounds(200, 100, 200, 25);
		pnlInputs.add(txtTestAllPath);
		
		JButton btnTestAll = new JButton("Test All");
		btnTestAll.setBounds(450, 100, 100, 25);
		pnlInputs.add(btnTestAll);
		btnTestAll.addActionListener
		(
				new ActionListener() {
					public void actionPerformed(ActionEvent e) 
					{
						if(isForAllDimensions)
						{
							runForAllDimensions();
						}
						else
						{
							onTestAllClicked();
						}
					}
				}
		);
		
		txtDebug.setBackground(DEBUG_BACKGROUND_COLOR);
		txtDebug.setBounds(MARGIN_FRAME, DEBUG_PANEL_YPOS, DEBUG_PANEL_WIDTH, DEBUG_PANEL_HEIGHT);
		txtDebug.setFont(new Font("Courier New", Font.PLAIN, 11));
		JScrollPane jspDebug = new JScrollPane(txtDebug);
		jspDebug.setBounds(MARGIN_FRAME, DEBUG_PANEL_YPOS, DEBUG_PANEL_WIDTH-15, DEBUG_PANEL_HEIGHT-15);
		mainContentPane.add(jspDebug);
		
		writeDebug("1. To train with specific dimension");
		writeDebug("    -> Enter dimension [1 to 100]");
		writeDebug("    -> Write path of folder where training images located");
		writeDebug("    -> Press Train button");
		writeDebug("2. To test single image after training");
		writeDebug("    -> Write path of image file");
		writeDebug("    -> Press Test button");
		writeDebug("3. To test images of a folder after training");
		writeDebug("    -> Write path of folder where test images located");
		writeDebug("    -> Press Test All button");
		writeDebug("3. To train all dimensions and test all images");
		writeDebug("    -> Enable 'For all dimensions' checkbox");
		writeDebug("    -> Write path of folder where training images located");
		writeDebug("    -> Write path of folder where test images located");
		writeDebug("    -> Press Test All button");
		
		txtTrainingImagePath.setText("train_images");
		txtTestImagePath.setText("test_images/s01/1.bmp");
		txtTestAllPath.setText("test_images");
	}

	private void runForAllDimensions() 
	{
		int [] failCounts = new int[100];
		int [] accuracies = new int[100];
		
		isDebugEnabled = false; // close debugs for faster results
		for(int i = 0; i < 100; i++)
		{
			spinFeatureDimension.setValue(i + 1);
			onTrainingClicked();
			failCounts[i] = onTestAllClicked();
			accuracies[i] = (int)(100 * (1 - ((double)failCounts[i] / (double)200)));
		}
		drawGraphichs(failCounts, accuracies);
		isDebugEnabled = true; // enable debugs for faster results
		
		String strAccuracies = "Accuracies : ";
		String strFailures = "Failures : ";
		for(int i = 0; i < 100; i++)
		{
			strAccuracies += String.format("%3d, ", accuracies[i]);
			strFailures += String.format("%3d, ", failCounts[i]);
		}
		writeDebug(strAccuracies);
		writeDebug(strFailures);
	}
	
	private void onTrainingClicked() 
	{
		// TODO Auto-generated method stub
		String strPath = txtTrainingImagePath.getText();
		
		java.io.File rootDir = new File(strPath);
		trainer = new AgentPCA((int)spinFeatureDimension.getValue());

		for (final File dirTrainset : rootDir.listFiles()) 
		{
	        if (dirTrainset.isDirectory()) 
	        {
	        	writeDebug(dirTrainset.getName() + " ... adding to training list");
	            for (final File fileImage : dirTrainset.listFiles()) 
	            {
	            	if (!fileImage.isDirectory()) 
	            	{
	            		if(!fileImage.getName().endsWith("bmp"))
	            		{
	            			continue;
	            		}
	            		trainer.addTrainingData(LabelMatrixPair.fromFile(fileImage));
	            		writeDebug("\t" + fileImage.getName() + " ... added to training list");
	            	}
	            }
	        } 
	    }
		trainer.train();
    	writeDebug("Training finished");
	}

	private void onTestClicked() 
	{
		LabelMatrixPair textMatrixPair= LabelMatrixPair.fromFile(new File(txtTestImagePath.getText()));
		LabelMatrixPair foundMatrixPair = trainer.findNearestNeighbor(textMatrixPair.getMatrix());

    	writeDebug("Search item : " + textMatrixPair.getLabel() + "; Found item : " + foundMatrixPair.getLabel());
	}

	private int onTestAllClicked() 
	{
		String strPath = txtTestAllPath.getText();
		
		java.io.File rootDir = new File(strPath);

		int totalCounter = 0;
		int failedCounter = 0;
		for (final File dirTestItems : rootDir.listFiles()) 
		{
	        if (dirTestItems.isDirectory()) 
	        {
	            for (final File fileTestImage : dirTestItems.listFiles()) 
	            {
	            	if (!fileTestImage.isDirectory() && fileTestImage.getName().endsWith("bmp")) 
	            	{
	            		LabelMatrixPair textMatrixPair= LabelMatrixPair.fromFile(fileTestImage);
	            		LabelMatrixPair foundMatrixPair = trainer.findNearestNeighbor(textMatrixPair.getMatrix());
	            		
	            		writeDebug("Search item : " + textMatrixPair.getLabel() + "; Found item : " + foundMatrixPair.getLabel());
	            		totalCounter++;
	            		if(!textMatrixPair.getLabel().equals(foundMatrixPair.getLabel()))
	            		{
	            			failedCounter++;
	            		}
	            	}
	            }
	        } 
	    }
    	writeDebug("Total : " + totalCounter + "; Failed : " + failedCounter);
    	
    	return failedCounter;
	}
	
	public void drawGraphichs(int [] failCounts, int [] accuracy)
	{
		int [] xValues = new int[100];
		for(int i = 0; i < 100; i++)
		{
			xValues[i] = i + 1;
		}

		Graph2D graphPanelAccuracy = new Graph2D("Accuracies");

		graphPanelAccuracy.setPoints(xValues, accuracy);
		graphPanelAccuracy.setLocation(this.getX() + 150, this.getY() + 150);
		graphPanelAccuracy.setSize(300, 300);
		graphPanelAccuracy.setVisible(true);
		/*
		Graph2D graphPanelFails = new Graph2D("Fail Counts");
		
		graphPanelFails.setPoints(xValues, failCounts);
		graphPanelFails.setLocation(this.getX() + 150, this.getY() + 150);
		graphPanelFails.setSize(300, 300);
		graphPanelFails.setVisible(true);*/
	}
	public void clearDebug() 
	{
		txtDebug.setText("");
	}
	public void writeDebug(String message) 
	{
		if(isDebugEnabled)
		{
			txtDebug.setText(txtDebug.getText() + "\n" + message );
		}
	}
	
	public static void main(String[] args) 
	{
		EventQueue.invokeLater(new Runnable() 
		{
			public void run() 
			{
				try 
				{					
					MainFrame frame = new MainFrame();
					frame.setVisible(true);
				} 
				catch (Exception e) 
				{
					e.printStackTrace();
				}
			}
		});

	}

}
