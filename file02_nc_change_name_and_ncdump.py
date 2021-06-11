# coding=UTF-8
import os
import string

def main(path):
    filenames = os.listdir('%s/HALO_DB/' %path)
    dates = ['20180308', '20180310', '20180312', '20180317', '20180319',\
             '20180320', '20180322', '20180324', '20180326', '20180328',\
             '20180330', '20180403', '20180404', '20180407', '20180409']
    for d in dates:
        if not os.path.exists('preprocess/%s' %d):
            os.popen('mkdir preprocess/%s' %d)
    for fullname in filenames:
        if all(['adlr' in fullname, 'nc' in fullname]):
            print('{}:'.format(fullname))
            name, _ = fullname.split('.')
            if '20180308' in name or '20180403a12' in name:
                _, dataset, _, obj, date_2 = name.split('_')
            else:
                _, dataset, _, obj, date_2, _ = name.split('_')
            date = date_2.replace('a12','').replace('a','')
            print(date)
            storage_name = '_'.join([obj, date, dataset]) 
            if not os.path.exists('%s/preprocess/%s/%s.nc' %(path, date, storage_name)):
                cmd1 = 'cp %s/HALO_DB/%s.nc %s/preprocess/%s/%s.nc' %(path, name, path, date, storage_name)
                os.popen(cmd1)
                cmd2 = 'ncdump -h %s/HALO_DB/%s.nc >& %s/preprocess/%s/%s_info.txt' %(path, name, path, date, storage_name)
                os.popen(cmd2)
                print('Finished!')
            else:
                print('File exists!')
        else:
            print('Skip '+fullname)

if __name__ == "__main__":
    main(path='.')