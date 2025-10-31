void visu()
{

	int run_number;
	int event_number;
	int solution;
	
	cout << "Choose run (-1 to show legend): ";
	cin  >> run_number; 
	
	if(run_number == -1)
	{
		cout << endl << "3D model contains:" << endl;
		cout << "	Calorimeter hits: low treshold hits (yellow)" << endl;
		cout << "	                  high treshold hits (red)" << endl;
		cout << "	Tracker hits: hits with usable radii (red)" << endl;
		cout << "	              hits with missing or wrong radii (yellow)" << endl;
		cout << "	              hits with missing height are set to z = 0.0" << endl;
		cout << "	Bi calibration sources: light blue (positions not clear at the moment)" << endl;
		cout << "	Reconstrcution (if calculated): one line per side of detector (blue)" << endl << endl;
	}
	else
	{
		cout << "Choose event to visualize (-1 to exit): ";
		cin  >> event_number; 
		
		cout << "Choose solution to visualize (-1 to exit): ";
		cin  >> solution; 
		
		
		if(event_number != -1)
		{
			TCanvas* c = new TCanvas("3D_demonstrator", "", 1920, 1080);
						
			TGLViewer* view = (TGLViewer*)gPad->GetViewer3D();
			TGeoManager* manager = new TGeoManager();
			
			TFile* file = new TFile(Form("./Events_visu/Run-%d_event-%d_solution-%d_3D.root", run_number, event_number, solution));
			TGeoVolume* geo;
			file->GetObject("demonstrator", geo);
			
			manager->SetTopVolume(geo);
			manager->CloseGeometry(); 
			geo->Draw();
			
			int i = 0;
			TPolyLine3D *track;
			while(file->GetListOfKeys()->Contains(Form("track-%d", i)))
			{
				file->GetObject(Form("track-%d", i), track);
				track->Draw("same");
				++i;
			}
			
			i = 0;
			TPolyLine3D *trajectory;
			while(file->GetListOfKeys()->Contains(Form("trajectory-%d", i)))
			{
				file->GetObject(Form("trajectory-%d", i), trajectory);
				trajectory->Draw("same");
				++i;
			}
			
			i = 0;
			TPolyMarker3D* trajectory_points;
			while(file->GetListOfKeys()->Contains(Form("trajectory_points-%d", i)))
			{
				file->GetObject(Form("trajectory_points-%d", i), trajectory_points);
				trajectory_points->Draw("same");
				++i;
			}
			
			TPolyMarker3D* association_points;
			file->GetObject( "association_points", association_points );
			association_points->Draw("same");
			
			TLatex* title = new TLatex(-0.9, 0.9, Form("Run %d | Event %d | Solution %d", run_number, event_number, solution));
			title->SetTextSize(0.03);
			title->Draw("same");
			
		}
	}
}


