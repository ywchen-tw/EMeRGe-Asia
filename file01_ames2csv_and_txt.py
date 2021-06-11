# coding=UTF-8
from __future__ import print_function
import nastools
import csv
import numpy as np
import os


def main():
    # get all original files
    # EMeRGe-Asia's data are downloaded from https://halo-db.pa.op.dlr.de/mission/97
    filenames = os.listdir('HALO_DB/')
    # make separate directories for each date
    dates = ['20180308','20180310','20180312','20180317','20180319','20180320','20180322','20180324',\
             '20180326','20180328','20180330','20180403','20180404','20180407','20180409']
    for date in dates:
        if not os.path.exists('preprocess/'):
            cmd = 'mkdir preprocess'
            os.popen(cmd)
        if not os.path.exists('preprocess/{}'.format(date)):
            cmd = 'mkdir preprocess/{}'.format(date)
            os.popen(cmd)


    for fullname in filenames:
        if 'adlr' not in fullname and '.ipynb' not in fullname and '.txt' not in fullname:    # exclude adlr files
            print(fullname+':') 
            name, _ = fullname.split('.')
            # parse the filename into separate elements
            if 'PANGC' in name or 'FAIROCI' in name or 'SO2' in name:
                _, dataset, _, _, _, date, instrument, obj, ver = name.split('_')
            elif 'CH4_CO2' in name:
                _, dataset, _, _, _, date, obj1, obj2, ver = name.split('_')
                obj = obj1 + obj2
                instrument = 'PICARRO'
            elif any(['AMETYST' in name]):
                _, dataset, _, _, _, _, date, instrument, obj, ver = name.split('_')
                date = date[:-1]
            elif '10VOCs' in name:
                #HALO-DB_dataset6052_release4_EMAS16_20180409_HKMS_10VOCs_V05.ames
                _, dataset, _, _, date, instrument, obj, ver = name.split('_')
            else:
                _, dataset, _, _, date, instrument, obj, ver = name.split('_')
            storage_name = '_'.join([obj, date, instrument, dataset]) 
            if not os.path.exists('preprocess/{}/{}.csv'.format(date, storage_name)):
                h = nastools.Naspy('HALO_DB/'+name+".ames")
                fp = open('preprocess/{}/{}.csv'.format(date, storage_name), "w")
                fi = open('preprocess/{}/{}_info.txt'.format(date, storage_name), "w")
                header_list= []
                for key, val in h.header.__dict__.iteritems():
                    fi.write('%s: %s\n' % (key, val)) 
                indep = h.header.INDEPENDENT_VARIABLE
                header_list.append(indep['NAME']+'('+indep['DESC']+')')
                dep = h.header.DEPENDENT_VARIABLE
                for d in dep:
                    header_list.append(d['NAME']+'('+d['DESC']+')')
                # Get variable names
                writer = csv.writer(fp)
                # print h.get_column_names
                writer.writerow(h.get_column_names())
                writer.writerow(header_list)
                # Make a masked array with missing values masked
                arr = h.make_numpy(masked=True)
                for row in arr:
                    row = list(row)
                    writer.writerow(row)


                # Make a pandas.DataFrame
                # df = h.make_DataFrame()

                # Make a pandas.DataFrame with integer index (not datetime)
                # df = h.make_DataFrame(datetime_asindex=False)
                # rename 'DATETIME' field to 'DT'
                # df = df.rename(columns={'DATETIME':'DT'})

                #fp.writelines(lines)
                print('finish {}!'.format(name))
                fp.close()
                fi.close()
            else:
                print('File exists!')
        else:
            print('Skip {}.'.format(fullname))

if __name__ == "__main__":
    main()