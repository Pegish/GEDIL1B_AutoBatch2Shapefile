# GEDIL1B_AutoBatch2Shapefile


  # Import Libraries
    import os
    import h5py
    import numpy as np
    import pandas as pd
    import geopandas as gp
    import shapefile
    from shapely.geometry import Point
    from shapely.geometry import box
    from shapely.geometry import MultiPolygon
    import geoviews as gv
    from geoviews import opts, tile_sources as gvts
    import holoviews as hv
    gv.extension('bokeh', 'matplotlib')
    import shapely
    import warnings
    from shapely.errors import ShapelyDeprecationWarning
    warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning) 
  
  # Verify the change and directory
    current_directory = os.getcwd()
    print("Current Directory:", current_directory)
    new_directory = "D:/Python/NewFilter-GEDIL1BGS2020"
    os.chdir(new_directory)
  
    inDir = os.getcwd()
    print("Updated Directory:", inDir)
    output_folder = os.path.join(new_directory, "outputpolygonForestwvNoSEM")
    os.makedirs(output_folder, exist_ok=True)
    print("Output Folder:", output_folder)
    
  # Specify the path to your shapefile of study area (ROI)
    shapefile_path = 'D:/Python/NewFilter-GEDIL1BGS2020/forestproject.shp'

  # Step 1: Read the shapefile using geopandas
    forest = gp.read_file(shapefile_path)

  # Step 2: Initialize lists to store lengths
    lengths_variable_selection = []
    lengths_after_filtering = []
    lengths_after_clipping = []
    
  # Step 3: List GEDI L2A .h5 files in the inDir
    gediFiles = [g for g in os.listdir(inDir) if g.startswith('processed_GEDI01_B') and g.endswith('.h5')]  
    
  # Step 4: Explore metadata and beam names
      for L1B_file in gediFiles:
          L1B_path = os.path.join(inDir, L1B_file)
          #print(f"Processing file: {L2A_path}")
      
          with h5py.File(L1B_path, 'r') as gediL1B:
             
              print("Keys in the file:", list(gediL1B.keys()))
              print("METADATA:", list(gediL1B['METADATA']))
              for g in gediL1B['METADATA']['DatasetIdentification'].attrs:
                  print(g)
              print("Purpose:", gediL1B['METADATA']['DatasetIdentification'].attrs['purpose'])
      
              beamNames = [g for g in gediL1B.keys() if g.startswith('BEAM')]
              #print("Beam Names:", beamNames)

   # Step 5: Initialize lists to store data
              gediL1B_objs = []
              gediL1B.visit(gediL1B_objs.append)
              gediSDS = [o for o in gediL1B_objs if isinstance(gediL1B[o], h5py.Dataset)]
              #gediSDS = [key for key in gediL1B.keys() if isinstance(gediL1B[key], h5py.Dataset)]
              
              shotNum, dem, demsrtm, zElevation, zHigh, zLat, zLon, srf, degrade, SolarElev, Solarazi, beamI, rxcount, rxstart = ([] for _ in range(14))
  # Step 6: Read data from each beam and append to lists
              for beamName in beamNames:
                 [shotNum.append(str(h)) for h in gediL1B[[g for g in gediSDS if g.endswith('/shot_number') and beamName in g][0]][()]]
                 [dem.append(h) for h in gediL1B[[g for g in gediSDS if g.endswith('/digital_elevation_model') and beamName in g][0]][()]]
                 [demsrtm.append(h) for h in gediL1B[[g for g in gediSDS if g.endswith('/digital_elevation_model_srtm') and beamName in g][0]][()]]
                 [zElevation.append(h) for h in gediL1B[[g for g in gediSDS if g.endswith('/elevation_bin0') and beamName in g][0]][()]]  
                 [zHigh.append(h) for h in gediL1B[[g for g in gediSDS if g.endswith('/elevation_lastbin') and beamName in g][0]][()]]  
                 [zLat.append(h) for h in gediL1B[[g for g in gediSDS if g.endswith('/latitude_bin0') and beamName in g][0]][()]]  
                 [zLon.append(h) for h in gediL1B[[g for g in gediSDS if g.endswith('/longitude_bin0') and beamName in g][0]][()]]  
                 [srf.append(h) for h in gediL1B[[g for g in gediSDS if g.endswith('/stale_return_flag') and beamName in g][0]][()]]  
                 [degrade.append(h) for h in gediL1B[[g for g in gediSDS if g.endswith('/degrade') and beamName in g][0]][()]]  
                 [SolarElev.append(h) for h in gediL1B[[g for g in gediSDS if g.endswith('/solar_elevation') and beamName in g][0]][()]]
                 [Solarazi.append(h) for h in gediL1B[[g for g in gediSDS if g.endswith('/solar_azimuth') and beamName in g][0]][()]]
                 [rxcount.append(h) for h in gediL1B[[g for g in gediSDS if g.endswith('/rx_sample_count') and beamName in g][0]][()]]
                 [rxstart.append(h) for h in gediL1B[[g for g in gediSDS if g.endswith('/rx_sample_start_index') and beamName in g][0]][()]]
                 [beamI.append(h) for h in [beamName] * len(gediL1B[[g for g in gediSDS if g.endswith('/shot_number') and beamName in g][0]][()])]
                 #[rxwv.append(h) for h in gediL1B[[g for g in gediSDS if g.endswith('/rxwaveform') and beamName in g][0]][()]]
                 #[txwv.append(h) for h in gediL1B[[g for g in gediSDS if g.endswith('/txwaveform') and beamName in g][0]][()]]
              
   # Step 7: Store the data in a DataFrame
                dataframe = pd.DataFrame({
                    'Shot Number': shotNum,
                    'Beam': beamI,
                    'Latitude': zLat,
                    'Longitude': zLon,
                    'Digital Elevation Model': dem,
                    'SRTM DEM': demsrtm,
                    'Ground Elevation (m)': zElevation,
                    'Canopy Elevation (m)': zHigh,
                    'Stale Return Flag': srf,
                    'Degrade Flag': degrade,
                    'Solar Elevation': SolarElev,
                    'Solar Azimuth': Solarazi,
                    'rx Sample Count': rxcount,
                    'rx Start Index': rxstart,
                    #'rxWaveform': rxwv,
                    #'txWaveform': txwv
                })
                lengths_variable_selection.append(len(dataframe))
                # Read RH data from HDF5 file and append to the list
                del shotNum, dem, demsrtm, zElevation, zHigh, zLat, zLon, srf, degrade, SolarElev, Solarazi, beamI
                print(len(dataframe))
                
   # Step 8: Import GeoJSON as GeoDataFrame
                ROI = gp.GeoDataFrame.from_file('ROI.geojson') 
                ROI ['geometry'][0]  # Plot GeoDataFrame      
                ROI.envelope[0].bounds
                minLon, minLat, maxLon, maxLat = ROI.envelope[0].bounds  # Define the min/max lat/lon from the bounds of Redwood NP    
                dataframe = dataframe.where(dataframe['Latitude'] > minLat)
                dataframe = dataframe.where(dataframe['Latitude'] < maxLat)
                dataframe = dataframe.where(dataframe['Longitude'] > minLon)
                dataframe = dataframe.where(dataframe['Longitude'] < maxLon)
                dataframe = dataframe.dropna()
                print(len(dataframe))
                print("dataframe after variable selection:", len(dataframe))
        
   # Step 9: Apply filters to remove specific variables
                dataframe = dataframe.where(dataframe['Stale Return Flag'].ne(1))
                dataframe = dataframe.where(dataframe['Degrade Flag'] < 1) 
                #dataframe = dataframe.where((dataframe['Digital Elevation Model']-dataframe['Ground Elevation (m)']) < abs(2))
                #dataframe = dataframe.where(dataframe['Solar Elevation'] < 0)
                dataframe = dataframe.dropna()
                lengths_after_filtering.append(len(dataframe))
                print(len(dataframe))
                #print(dataframe.dtypes)
                print("dataframe after filtering:", len(dataframe))
        
                waves = []  # List to store each list of waveform values
                for _, shot in dataframe.iterrows():
                    starts = int(shot['rx Start Index'])
                    wave = gediL1B[f"{shot['Beam']}/rxwaveform"][starts: starts + int(shot['rx Sample Count'])]
                    
                    # In order to export as GeoJSON, the type needs to be converted from a numpy array to a single comma separated string
                    waves.append(','.join([str(q) for q in wave]))
            
                # Add to data frame
                dataframe['rx Waveform'] = waves
        
   # Step 10: Convert the DataFrame to GeoDataFrame
                geometry = [Point(xy) for xy in zip(dataframe['Longitude'], dataframe['Latitude'])]
                crs = {'init': 'epsg:4326'}
                gdf = gp.GeoDataFrame(dataframe, crs=crs, geometry=geometry)
        
                # Create circles around each point with a radius of 12.5 meters
                radius = 0.0001126126  # Radius in meters
                gdf['geometry'] = gdf['geometry'].buffer(radius)
        
                gdf_forest = gp.sjoin(gdf, forest, how="inner", op="intersects")
                lengths_after_clipping.append(len(gdf_forest))
  # Step 11: Save the GeoDataFrame with polygons within the forest as shapefile
                out_shapefile = os.path.join(output_folder, L1B_file.replace('.h5', '_forestwvNoDEM.shp'))
                gdf_forest.to_file(out_shapefile)
        
                # Save the GeoDataFrame with polygons within the forest as GeoJSON
                out_geojson = os.path.join(output_folder, L1B_file.replace('.h5', '_forestwvNoDEM.json'))
                gdf_forest.to_file(out_geojson, driver='GeoJSON')
                # Clean up
                del dataframe, gdf, waves


  # Step 12: After processing all files, create a dataframe from the lists
    lengths_df = pd.DataFrame({
        'Length after variable selection': lengths_variable_selection,
        'Length after filtering': lengths_after_filtering, 
        'Length after clipping' : lengths_after_clipping
    })

  # Step 13: Save the dataframe to a CSV file
    lengths_df.to_csv('dataframe_lengthsForestwvNoDEM.csv', index=False)

