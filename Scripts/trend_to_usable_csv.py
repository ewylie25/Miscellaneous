#! /usr/bin/python
"""
Script to convert unaligned csv exports of time seried data to a csv with a 
single date-time variable. Should work with arbitrary number of points.
March 2015
"""

from datetime import datetime as DM
import time as Time
import argparse

#Some input CSV formatting values for if this changes
count_heading_rows = 2
row_with_point_names = 1  # indexed from 0
count_columns_per_point = 4
format_point_data = {
    'date_column_index': 0, # indexed from 0
    'date_format': '%m/%d/%Y',
    'time_column_index': 1, # indexed from 0
    'time_format': '%I:%M:%S %p',
    'value_column_index': -1 # indexed from 0
}
empty_data_point = [''] * count_columns_per_point


def chunks(list_of_row, size):
    """ Yield successive n-sized chunks from provided row data."""
    for index in xrange(0, len(list_of_row), size):
        yield list_of_row[index : index + size]


def convert_to_datetime(old_data):
    """Convert strings of date and time to datetime objects in data structure.
        Must be modified if format changes."""
    data = []
    times = set()
    d = format_point_data['date_column_index']
    t = format_point_data['time_column_index']
    v = format_point_data['value_column_index']
    dt_string = " ".join([format_point_data['date_format'],format_point_data['time_format']])
    for row in old_data:
        new_row = [[DM.strptime(" ".join([point[d],point[t]]), dt_string), point[v]] if not any([point == empty_data_point, point == ['']]) else ['',''] for point in row]
        times.update(tuple([data_point[0] for data_point in new_row if data_point[0] != '']))
        data.append(new_row)
    return data, times


if __name__ == '__main__':
    # Dealing with arguments
    parser = argparse.ArgumentParser(description="Script to convert unaligned csv exports of time seried data", epilog="Contact Liz Wylie, ewylie25@gmail.com for help.")
    parser.add_argument('-i', '--input_file', type=str, required=True,help="input file fromatted by eDNA Trend")
    parser.add_argument('-o', '--output_file', type=str,help="filename to write output csv, otherwise defaults to './output_%Y%m%d_%H%M.csv' ")
    args = parser.parse_args()
    f_input = args.input_file
    f_output = args.output_file if args.output_file else './output_{0}.csv'.format(Time.strftime('%Y%m%d_%H%M'))    

    # read in data
    with open(f_input) as f:
        raw_text = f.readlines()

    
    # pull out point names and structure data into groups of point data
    name_points = [name for name in raw_text[row_with_point_names].strip('\n').split(",") if name!='']
    data = [list(chunks(row.strip('\n').split(","), count_columns_per_point)) for row in raw_text[count_heading_rows:]]

    # clean up data - convert to datetimes, get set of all datetimes
    data, times = convert_to_datetime(data)

    # assign data to  point names, gives one tuple per point
    formatted_data= zip(name_points, *data)

    # sort times
    all_times = sorted(list(times))

    #convert data columns to dictionary of {time:value}
    data_as_dict = []
    for column in formatted_data:
        time_value_dict = {data_point[0]:data_point[1] for data_point in column[1:]}
        new_column = tuple([column[0], time_value_dict])
        data_as_dict.append(new_column)

    #create lines in csv
    f_to_write = []
    f_to_write.append(",".join(["Datetime"]+[i[0] for i in data_as_dict]))
    for time in all_times:
        line = []
        line.append(time.strftime('%m/%d/%Y %I:%M:%S %p'))
        for i in data_as_dict:
            line.append(i[1].get(time, ""))
        f_to_write.append(",".join(line))

    #write to file
    with open(f_output, 'w') as f:
        f.write("\n".join(f_to_write))




