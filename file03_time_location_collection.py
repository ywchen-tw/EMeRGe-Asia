import csv
import HALO_processing as HALO

def date_processing_west_coast(date):
    print('Processing now (west coast):', date)
    if HALO.get_adlr_data(date) != None:
        # read nc file
        t_local_f, lat, lon, height, u, v, w, p_H, static_P, dynamic_P, TAT, TD, THETA_aircraft, THETA_V, TV, TS,\
        RELHUM, MIXRATIOV_H2O, MIXRATIOM_H2O = HALO.get_adlr_data(date)         # t_local for flight of nc file
        t_range = [0] if HALO.in_western_Taiwan_region(lat[0], lon[0]) else []
        for i in range(len(lat)-1):
            flag1 = True if HALO.in_western_Taiwan_region(lat[i], lon[i]) else False
            flag2 = True if HALO.in_western_Taiwan_region(lat[i+1], lon[i+1]) else False
            if flag1 != flag2:
                t_range.append(i)
        if HALO.in_western_Taiwan_region(lat[-1], lon[-1]):
            t_range.append(len(lat)-1)
        fp2 = open('preprocess/others/t_segment_western_Taiwan.csv', 'a')
        writer2 = csv.writer(fp2)
        t_range_f = t_range
        for t in t_range_f:
            content = [date, t_local_f[t]]
            writer2.writerow(content)
        fp2.close()
        if len(t_range_f) > 0:
            second_to_average = 60  # YCC lab (west coast)
            nc_t_interval = 10*second_to_average+1
            u_ma = HALO.moving_average(u, n=nc_t_interval)
            v_ma = HALO.moving_average(v, n=nc_t_interval)
            w_ma = HALO.moving_average(w, n=nc_t_interval)
            count = 0
            while count < len(t_range_f)-1:
                t_run_f = HALO.t_interval(t_range_f, count)
                if HALO.in_western_Taiwan_region(lat[t_run_f[2]], lon[t_run_f[2]]):
                    d = str(date)
                    mmyy = d.replace('2018', '')
                    fp = open('preprocess/others/loc_time_info_western_Taiwan.csv', 'a')
                    writer = csv.writer(fp)
                    t = range(t_run_f[0], t_run_f[-1])
                    for i in t:
                        conctent = [int(date), t_local_f[i], lat[i], lon[i], height[i], round(height[i]/100)*100,\
                                    u[i], v[i], w[i], p_H[i], static_P[i], dynamic_P[i], TAT[i], TD[i],\
                                    THETA_aircraft[i], THETA_V[i], TV[i], TS[i], RELHUM[i], MIXRATIOV_H2O[i], MIXRATIOM_H2O[i]]
                        writer.writerow(conctent)
                    fp.close()
                count += 1
    else:
        print('File {} does not exist!'.format(date))

def date_processing_all_area(date):
    print('Processing now (all area):', date)
    if HALO.get_adlr_data(date) != None:
        # read nc file
        t_local_f, lat, lon, height, u, v, w, p_H, static_P, dynamic_P, TAT, TD, THETA_aircraft, THETA_V, TV, TS,\
        RELHUM, MIXRATIOV_H2O, MIXRATIOM_H2O = HALO.get_adlr_data(date)         # t_local for flight of nc file
        t_range = [0, len(lat)-1]
        fp2 = open('preprocess/others/t_segment_all_area.csv', 'a')
        writer2 = csv.writer(fp2)
        t_range_f = t_range
        for t in t_range_f:
            content = [date, t_local_f[t]]
            writer2.writerow(content)
        fp2.close()
        second_to_average = 15  # YCC lab (all area)
        nc_t_interval = 10*second_to_average+1
        u_ma = HALO.moving_average(u, n=nc_t_interval)
        v_ma = HALO.moving_average(v, n=nc_t_interval)
        w_ma = HALO.moving_average(w, n=nc_t_interval)
        count = 0
        while count < len(t_range_f)-1:
            t_run_f = HALO.t_interval(t_range_f, count)
            d = str(date)
            mmyy = d.replace('2018', '')
            fp = open('preprocess/others/loc_time_info_all_area.csv', 'a')
            writer = csv.writer(fp)
            t = range(t_run_f[0], t_run_f[-1])
            for i in t:
                conctent = [int(date), t_local_f[i], lat[i], lon[i], height[i], round(height[i]/100)*100,\
                            u[i], v[i], w[i], p_H[i], static_P[i], dynamic_P[i], TAT[i], TD[i],\
                            THETA_aircraft[i], THETA_V[i], TV[i], TS[i], RELHUM[i], MIXRATIOV_H2O[i], MIXRATIOM_H2O[i]]
                writer.writerow(conctent)
            fp.close()
            count += 1
    else:
        print('File {} does not exist!'.format(date))

def main():
    date_set = ['20180308', '20180310', '20180312', '20180317', '20180319', '20180320', '20180322', '20180324', \
                '20180326', '20180328', '20180330', '20180403', '20180404', '20180407', '20180409']
    all_loc_header = ['date', 'local_time_sec', 'lat', 'lon', 'H', 'H_round_to_hundred', 'u_moving_average', 'v_moving_average',\
              'w_moving_average', 'pressure_H', 'static_P', 'dynamic_P', 'Total_Air_Temp', 'TD',\
              'THETA_aircraft', 'THETA_V', 'TV', 'TS','RELHUM', 'MIXRATIOV_H2O', 'MIXRATIOM_H2O']
    if not os.path.exists('preprocess/others'):
        cmd = 'mkdir preprocess/others'
        os.popen(cmd)
    # make file
    loc_time_info_files = ['preprocess/others/loc_time_info_western_Taiwan.csv',
                           'preprocess/others/loc_time_info_all_area.csv']
    t_segment_files = ['preprocess/others/t_segment_western_Taiwan.csv',
                       'preprocess/others/t_segment_all_area.csv']
    for file in loc_time_info_files:
        fp = open(file, 'w')
        writer = csv.writer(fp)
        writer.writerow(all_loc_header)
        fp.close()
    for file in t_segment_files:
        fp2 = open(file, 'w')
        writer2 = csv.writer(fp2)
        header2 = ['date', 'local_time_sec']
        writer2.writerow(header2)
        fp2.close()
    for d in date_set:
        date_processing_west_coast(d)
        date_processing_all_area(d)
        print('-'*10)
    print('All files are processed!!!')

if __name__ == "__main__":
    main()
