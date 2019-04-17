import pandas as pd
import numpy as np





def convert_bamberg_to_timeseries(
    dataframe, 
    time_line,
    pain_duration = 5,
    column_to_use_for_time = "Time_Relative_sf", 
    results_col_name_for_feature ="feature"
    ):
    """
    Convert the bamberg results to a timeseries

    Keyword arguments:
    dataframe -- contains the Bamberg data(default 0.0)
    time_line -- list with times onto which the bamberg data will be converted 
    pain_duration -- time in seconds after the start(default 5.0)
    column_to_use_for_time --
    results_col_name_for_feature --
    """

    key_time_rel = column_to_use_for_time
    key_feature = "Behavior"
    key_feature_type =  "Event_Type"
    key_feature_type_start = "State start"
    key_feature_type_stop = "State stop"
    key_feature_type_state = "State point"
    key_modifier = "Modifier_1"

    # Select the columns we are interested in
    selected_content = dataframe[[key_time_rel, key_feature, key_feature_type, key_modifier]]

    # The way the human FACS annotations are stored is in a way that the start and end of 
    # of an action unit is simply marked per line.
    stamps = []
    time_start, time_end, au, mod_start, mod_end = None, None, None, None, None
    

    # With i we walk over the rows of the dataframe
    for i in selected_content.index:
        # If a line contains a start of an AU we set the time and au ...
        if selected_content[key_feature_type][i] == key_feature_type_start and not time_start:
            time_start = selected_content[key_time_rel][i] 
            au = selected_content[key_feature][i]
            mod_start = selected_content[key_modifier][i]
            #print(i, selected_content["Event_Type"][i], selected_content["Behavior"][i])

            # ... and then we start to walk over the rest of the rows to find the end 
            # marker of this specific entry
            j = i
            while not time_end:
                #print(j)
                # Then if the end marker is found it is stored
                if (selected_content[key_feature_type][j] == key_feature_type_stop) and \
                        (selected_content[key_feature][j] == au):
                    time_end = selected_content[key_time_rel][j]
                    mod_end = selected_content[key_modifier][i]
                j+=1

            # If the end marker is found all is added to the stamps list
            if time_end:
                stamps.append([time_start, time_end, au, mod_start, mod_end])
                time_start, time_end, au, mod_start, mod_end = None, None, None, None, None
            # If the end is not found, we have an incomplete entry in the FACS file
            else:
                print("END NOT FOUND FOR: ", i, time_start, au)
                
        elif selected_content[key_feature_type][i] == key_feature_type_state:
            stamps.append([selected_content[key_time_rel][i], selected_content[key_time_rel][i]+pain_duration,\
                          selected_content[key_feature][i], np.nan, np.nan])
            
    
    start, end, features, modifier, modifier_end = zip(*stamps)

    # We now put our results in a dataframe
    timestamps = pd.DataFrame({"start": start, 
                           "end": end, 
                           results_col_name_for_feature: features,
                           "modifier": modifier})

    # Get the unique feature names
    unique_features = list(set(timestamps[results_col_name_for_feature]))
    # Make a new dataframe for the timeseries, using the features as columns. Set elements to zero
    timeseries = pd.DataFrame(data=0, index = np.arange(len(time_line)), columns=["time"]+unique_features)
    # Set the time column to the time_line given by the user
    timeseries["time"] = time_line

    # Now walk over all the detections in the timestamps dataframe
    for i, row in timestamps.iterrows():
        # Now set the elements in the time serie to the intensity at the time
        # of detection.
        timeseries[row[results_col_name_for_feature]].mask( 
            (timeseries["time"] >= row["start"]) & (timeseries["time"] <= row["end"]), 
            other=1 if np.isnan(row["modifier"]) else row["modifier"], 
            inplace=True )

    return timeseries
