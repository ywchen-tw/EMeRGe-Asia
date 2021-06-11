import csv
import time, datetime
from time import strftime, localtime
import numpy as np
import bisect
import HALO_processing as HALO
import glob, os

def t_ind_search(t_sec_obj_all, t):
    result = bisect.bisect_left(t_sec_obj_all, t)
    if result >= len(t_sec_obj_all):
        return len(t_sec_obj_all)-1
    else:
        return result

def t_start_end(t_sec_obj_all, t_now, t_formmer, t_latter):
    ind_now = t_ind_search(t_sec_obj_all, t_now)
    ind_now_f = t_ind_search(t_sec_obj_all, t_formmer)
    ind_now_l = t_ind_search(t_sec_obj_all, t_latter)
    return ind_now, ind_now_f, ind_now_l

def cal_nanmean(t_sec_obj_all, obj_all, t_now):
    t_formmer = t_now - 14
    t_latter = t_now + 1
    ind_now, ind_now_f, ind_now_l = t_start_end(t_sec_obj_all, t_now, t_formmer, t_latter)
    if all([abs(t_sec_obj_all[ind_now_f]-t_now) < 180., abs(t_sec_obj_all[ind_now_l]-t_now) < 180.]):
        return np.nanmean(np.array(obj_all[ind_now_f:ind_now_l])) 
    else:
        return np.nan

def main():
    dates = ['20180308', '20180310', '20180312', '20180317', '20180319', '20180320', '20180322', '20180324',\
             '20180326', '20180328', '20180330', '20180403', '20180404', '20180407', '20180409']
    date_nc, t_sec, lat, lon, H, H_round_to_hundred, u_ma, v_ma, w_ma,\
    p_H, static_P, dynamic_P, TAT, TD, THETA_aircraft, THETA_V, TV, TS,\
    RELHUM, MIXRATIOV_H2O, MIXRATIOM_H2O = HALO.adlr_data_csv('preprocess/others/loc_time_info_all_area.csv') 

    date_segment, t_seg = HALO.t_segment('preprocess/others/t_segment_all_area.csv') 
    t_sec_int = [round(t) for t in t_seg]
    t_ind = [bisect.bisect_left(t_sec,t) for t in t_sec_int]
    t_ind[-1] = t_ind[-1] - 1
    t_sec_SO2_all, SO2_all = ([] for _ in range(2))
    t_sec_CO_all, CO_all = ([] for _ in range(2))
    t_sec_NONOy_all, NO_all, NOy_all = ([] for _ in range(3))
    t_sec_O3_AMTEX_all, O3_AMTEX_all = ([] for _ in range(2))
    t_sec_O3_FAIROCI_all, O3_FAIROCI_all = ([] for _ in range(2))
    t_sec_CH4CO2_all, CH4_all, CO2_all = ([] for _ in range(3))
    t_sec_ORGNO3SO4NH4CHL_all, ALT_ORGNO3SO4NH4CHL_all, LAT_ORGNO3SO4NH4CHL_all, LON_ORGNO3SO4NH4CHL_all = ([] for _ in range(4))
    ORG_all, NO3_all, SO4_all, NH4_all, CHL_all = ([] for _ in range(5))
    ORG_PREC_all, NO3_PREC_all, SO4_PREC_all, NH4_PREC_all, CHL_PREC_all = ([] for _ in range(5))
    ORG_DL_all, NO3_DL_all, SO4_DL_all, NH4_DL_all, CHL_DL_all = ([] for _ in range(5))
    t_sec_PAN_all, LAT_PAN_all, LON_PAN_all, ALT_PAN_all, THETA_PAN_all, PAN_PPB_all = ([] for _ in range(6))
    t_sec_VOC_all = []
    FOR_all, FOR_PREC_all, FOR_LOD_all, MET_all, MET_PREC_all, MET_LOD_all = ([] for _ in range(6))
    ACN_all, ACN_PREC_all, ACA_all, ACA_PREC_all, ACE_all, ACE_PREC_all = ([] for _ in range(6))
    ISO_all, ISO_PREC_all = ([] for _ in range(2))
    BEN_all, BEN_PREC_all, TOL_all, TOL_PREC_all, XYL_all, XYL_PREC_all = ([] for _ in range(6))
    DMS_all, DMS_PREC_all, DMS_LOD_all, MVK_all, MVK_PREC_all, MEK_all, MEK_PREC_all = ([] for _ in range(7))
    t_sec_NO2_MINIDOAS_all, NO2_MINIDOAS_all, NO2_MINIDOAS_ERROR_all = ([] for _ in range(3))

    for d in dates:
        print('-'*50)
        print(d)
        if HALO.SO2_data_input(d) != None:
            t_sec_SO2, SO2 = HALO.SO2_data_input(d)
            print('SO2:', strftime("%m-%d %H:%M:%S", localtime(t_sec_SO2[0])))
            print('SO2:', strftime("%m-%d %H:%M:%S", localtime(t_sec_SO2[-1])))
            t_sec_SO2_all += list(t_sec_SO2)
            SO2_all += list(SO2)
            del t_sec_SO2, SO2
        if HALO.CO_data_input(d) != None:
            t_sec_CO, CO = HALO.CO_data_input(d)
            print('CO:', strftime("%m-%d %H:%M:%S", localtime(t_sec_CO[0])))
            print('CO:', strftime("%m-%d %H:%M:%S", localtime(t_sec_CO[-1])))
            t_sec_CO_all += list(t_sec_CO)
            CO_all += list(CO)
            del t_sec_CO, CO
        if HALO.NONOy_data_input(d) != None:
            t_sec_NONOy, NO, NOy = HALO.NONOy_data_input(d)
            print('NONOy:', strftime("%m-%d %H:%M:%S", localtime(t_sec_NONOy[0])))
            print('NONOy:', strftime("%m-%d %H:%M:%S", localtime(t_sec_NONOy[-1])))
            t_sec_NONOy_all += list(t_sec_NONOy)
            NO_all += list(NO)
            NOy_all += list(NOy)
            del t_sec_NONOy, NO, NOy
        if HALO.O3_AMTEX_data_input(d) != None:
            t_sec_O3_AMTEX, O3_AMTEX = HALO.O3_AMTEX_data_input(d)
            print('O3_AMTEX:', strftime("%m-%d %H:%M:%S", localtime(t_sec_O3_AMTEX[0])))
            print('O3_AMTEX:', strftime("%m-%d %H:%M:%S", localtime(t_sec_O3_AMTEX[-1])))
            t_sec_O3_AMTEX_all += list(t_sec_O3_AMTEX)
            O3_AMTEX_all += list(O3_AMTEX)
            del t_sec_O3_AMTEX, O3_AMTEX
        if HALO.O3_FAIROCI_data_input(d) != None:
            t_sec_O3_FAIROCI, O3_FAIROCI = HALO.O3_FAIROCI_data_input(d)
            print('O3_FAIROCI:', strftime("%m-%d %H:%M:%S", localtime(t_sec_O3_FAIROCI[0])))
            print('O3_FAIROCI:', strftime("%m-%d %H:%M:%S", localtime(t_sec_O3_FAIROCI[-1])))
            t_sec_O3_FAIROCI_all += list(t_sec_O3_FAIROCI)
            O3_FAIROCI_all += list(O3_FAIROCI)
            del t_sec_O3_FAIROCI, O3_FAIROCI
        if HALO.CH4CO2_data_input(d) != None:
            t_sec_CH4CO2, CH4, CO2 = HALO.CH4CO2_data_input(d)
            print('CH4CO2:', strftime("%m-%d %H:%M:%S", localtime(t_sec_CH4CO2[0])))
            print('CH4CO2:', strftime("%m-%d %H:%M:%S", localtime(t_sec_CH4CO2[-1])))
            t_sec_CH4CO2_all += list(t_sec_CH4CO2)
            CH4_all += list(CH4)
            CO2_all += list(CO2)
            del t_sec_CH4CO2, CH4, CO2
        if HALO.ORGNO3SO4NH4CHL_data_input(d) != None:
            t_sec_ORGNO3SO4NH4CHL, ALT_ORGNO3SO4NH4CHL, LAT_ORGNO3SO4NH4CHL, LON_ORGNO3SO4NH4CHL,\
                ORG, ORG_PREC, ORG_DL, NO3, NO3_PREC, NO3_DL, SO4, SO4_PREC, SO4_DL, NH4, NH4_PREC, NH4_DL,\
                CHL, CHL_PREC, CHL_DL = HALO.ORGNO3SO4NH4CHL_data_input(d)
            print('ORGNO3SO4NH4CHL:', strftime("%m-%d %H:%M:%S", localtime(t_sec_ORGNO3SO4NH4CHL[0])))
            print('ORGNO3SO4NH4CHL:', strftime("%m-%d %H:%M:%S", localtime(t_sec_ORGNO3SO4NH4CHL[-1])))
            t_sec_ORGNO3SO4NH4CHL_all += list(t_sec_ORGNO3SO4NH4CHL)
            ALT_ORGNO3SO4NH4CHL_all += list(ALT_ORGNO3SO4NH4CHL)
            LAT_ORGNO3SO4NH4CHL_all += list(LAT_ORGNO3SO4NH4CHL)
            LON_ORGNO3SO4NH4CHL_all += list(LON_ORGNO3SO4NH4CHL)
            ORG_all += list(ORG)
            NO3_all += list(NO3)
            SO4_all += list(SO4)
            NH4_all += list(NH4)
            CHL_all += list(CHL)
            ORG_PREC_all += list(ORG_PREC)
            NO3_PREC_all += list(NO3_PREC)
            SO4_PREC_all += list(SO4_PREC)
            NH4_PREC_all += list(NH4_PREC)
            CHL_PREC_all += list(CHL_PREC)
            ORG_DL_all += list(ORG_DL)
            NO3_DL_all += list(NO3_DL)
            SO4_DL_all += list(SO4_DL)
            NH4_DL_all += list(NH4_DL)
            CHL_DL_all += list(CHL_DL)
            del t_sec_ORGNO3SO4NH4CHL, ALT_ORGNO3SO4NH4CHL, LAT_ORGNO3SO4NH4CHL, LON_ORGNO3SO4NH4CHL, ORG, ORG_PREC, ORG_DL,\
                NO3, NO3_PREC, NO3_DL, SO4, SO4_PREC, SO4_DL, NH4, NH4_PREC, NH4_DL, CHL, CHL_PREC, CHL_DL 
        if HALO.PAN_data_input(d) != None:
            t_sec_PAN, LAT_PAN, LON_PAN, ALT_PAN, THETA_PAN, PAN_PPB = HALO.PAN_data_input(d)
            print('PAN:', strftime("%m-%d %H:%M:%S", localtime(t_sec_PAN[0])))
            print('PAN:', strftime("%m-%d %H:%M:%S", localtime(t_sec_PAN[-1])))
            t_sec_PAN_all += list(t_sec_PAN)
            LAT_PAN_all += list(LAT_PAN)
            LON_PAN_all += list(LON_PAN)
            ALT_PAN_all += list(ALT_PAN)
            THETA_PAN_all += list(THETA_PAN)
            PAN_PPB_all += list(PAN_PPB)
            del t_sec_PAN, LAT_PAN, LON_PAN, ALT_PAN, THETA_PAN, PAN_PPB
        if HALO.VOC_data_input(d) != None:
            t_sec_VOC, FOR, FOR_PREC, FOR_LOD, MET, MET_PREC, MET_LOD, ACN, ACN_PREC, ACA, ACA_PREC,\
                ACE, ACE_PREC, ISO, ISO_PREC, BEN, BEN_PREC, TOL, TOL_PREC, XYL, XYL_PREC,\
                DMS, DMS_PREC, DMS_LOD, MVK, MVK_PREC, MEK, MEK_PREC = HALO.VOC_data_input(d)
            print('VOC:', strftime("%m-%d %H:%M:%S", localtime(t_sec_VOC[0])))
            print('VOC:', strftime("%m-%d %H:%M:%S", localtime(t_sec_VOC[-1])))
            t_sec_VOC_all += list(t_sec_VOC)
            FOR_all += list(FOR)
            FOR_PREC_all += list(FOR_PREC)
            FOR_LOD_all += list(FOR_LOD)
            MET_all += list(MET)
            MET_PREC_all += list(MET_PREC)
            MET_LOD_all += list(MET_LOD)
            ACN_all += list(ACN)
            ACN_PREC_all += list(ACA_PREC)
            ACA_all += list(ACA)
            ACA_PREC_all += list(ACA_PREC)
            ACE_all += list(ACE)
            ACE_PREC_all += list(ACE_PREC)
            ISO_all += list(ISO)
            ISO_PREC_all += list(ISO_PREC)
            BEN_all += list(BEN)
            BEN_PREC_all += list(BEN_PREC)
            TOL_all += list(TOL)
            TOL_PREC_all += list(TOL_PREC)
            XYL_all += list(XYL)
            XYL_PREC_all += list(XYL_PREC)
            DMS_all += list(DMS)
            DMS_PREC_all += list(DMS_PREC)
            DMS_LOD_all += list(DMS_LOD)
            MVK_all += list(MVK)
            MVK_PREC_all += list(MVK_PREC)
            MEK_all += list(MEK)
            MEK_PREC_all += list(MEK_PREC)
            del t_sec_VOC, FOR, FOR_PREC, FOR_LOD, MET, MET_PREC, MET_LOD, ACN, ACN_PREC, ACA, ACA_PREC,\
                ACE, ACE_PREC, ISO, ISO_PREC, BEN, BEN_PREC, TOL, TOL_PREC, XYL, XYL_PREC,\
                DMS, DMS_PREC, DMS_LOD, MVK, MVK_PREC, MEK, MEK_PREC
        if HALO.NO2_MINIDOAS_data_input(d) != None:
            t_sec_NO2_MINIDOAS, NO2_MINIDOAS, NO2_MINIDOAS_ERROR = HALO.NO2_MINIDOAS_data_input(d)
            print('NO2_MINIDOAS:', strftime("%m-%d %H:%M:%S", localtime(t_sec_NO2_MINIDOAS[0])))
            print('NO2_MINIDOAS:', strftime("%m-%d %H:%M:%S", localtime(t_sec_NO2_MINIDOAS[-1])))
            t_sec_NO2_MINIDOAS_all += list(t_sec_NO2_MINIDOAS)
            NO2_MINIDOAS_all += list(NO2_MINIDOAS)
            NO2_MINIDOAS_ERROR_all += list(NO2_MINIDOAS_ERROR)
            del t_sec_NO2_MINIDOAS, NO2_MINIDOAS, NO2_MINIDOAS_ERROR

    fp = open('Analysis_dataset/All_area_dataset-15s.csv','w')
    writer = csv.writer(fp)
    header = ['date', 'Local time', 'UTC_time', 'lat', 'lon', 'height', 'height_in_hundred', 'u_ma', 'v_ma', 'w_ma',\
              'pressure_H', 'static_P', 'dynamic_P', 'Total_Air_Temp', 'TD', 'THETA_aircraft', 'THETA_V', 'TV', 'TS',\
              'RELHUM', 'MIXRATIOV_H2O', 'MIXRATIOM_H2O', 'SO2', 'CO', 'NO', 'NOy', 'O3_AMTEX','O3_FAIROCI', 'CH4', 'CO2',\
              'ORG', 'ORG_PREC', 'ORG_DL', 'NO3', 'NO3_PREC', 'NO3_DL', 'SO4', 'SO4_PREC', 'SO4_DL',\
              'NH4', 'NH4_PREC', 'NH4_DL', 'CHL', 'CHL_PREC', 'CHL_DL',\
              'LAT_PAN', 'LON_PAN', 'ALT_PAN', 'THETA_PAN', 'PAN_PPB', 'FOR', 'FOR_PREC', 'FOR_LOD',\
              'MET', 'MET_PREC', 'MET_LOD', 'ACN', 'ACN_PREC', 'ACA', 'ACA_PREC', 'ACE', 'ACE_PREC',\
              'ISO', 'ISO_PREC',  'BEN', 'BEN_PREC', 'TOL', 'TOL_PREC', 'XYL', 'XYL_PREC',\
              'DMS', 'DMS_PREC', 'DMS_LOD', 'MVK', 'MVK_PREC',  'MEK', 'MEK_PREC',\
              'NO2_MINIDOAS', 'NO2_MINIDOAS_ERROR',]
    writer.writerow(header)
    for i in range(len(t_sec_int)-1):
        for j in range(t_ind[i],t_ind[i+1]):
            t_now = t_sec[j]
            if t_sec[j]%15 == 0:
                t_formmer = t_now - 14
                t_latter = t_now + 1
                # SO2
                SO2_write = cal_nanmean(t_sec_SO2_all,SO2_all,t_now)
                # CO
                CO_write = cal_nanmean(t_sec_CO_all,CO_all,t_now)
                # NONOy
                NO_write = cal_nanmean(t_sec_NONOy_all,NO_all,t_now)
                NOy_write = cal_nanmean(t_sec_NONOy_all,NOy_all,t_now)
                # O3_AMTEX
                O3_AMTEX_write = cal_nanmean(t_sec_O3_AMTEX_all,O3_AMTEX_all,t_now)
                # O3_FAIROCI
                O3_FAIROCI_write = cal_nanmean(t_sec_O3_FAIROCI_all,O3_FAIROCI_all,t_now)
                # CH4CO2
                CH4_write = cal_nanmean(t_sec_CH4CO2_all,CH4_all,t_now)
                CO2_write = cal_nanmean(t_sec_CH4CO2_all,CO2_all,t_now)
                # ORGNO3SO4NH4CHL, remove negative value
                ORG_write = cal_nanmean(t_sec_ORGNO3SO4NH4CHL_all,ORG_all,t_now)
                ORG_PREC_write = cal_nanmean(t_sec_ORGNO3SO4NH4CHL_all,ORG_PREC_all,t_now)
                ORG_DL_write = cal_nanmean(t_sec_ORGNO3SO4NH4CHL_all,ORG_DL_all,t_now)
                NO3_write = cal_nanmean(t_sec_ORGNO3SO4NH4CHL_all,NO3_all,t_now)
                NO3_PREC_write = cal_nanmean(t_sec_ORGNO3SO4NH4CHL_all,NO3_PREC_all,t_now)
                NO3_DL_write = cal_nanmean(t_sec_ORGNO3SO4NH4CHL_all,NO3_DL_all,t_now)
                SO4_write = cal_nanmean(t_sec_ORGNO3SO4NH4CHL_all,SO4_all,t_now)
                SO4_PREC_write = cal_nanmean(t_sec_ORGNO3SO4NH4CHL_all,SO4_PREC_all,t_now)
                SO4_DL_write = cal_nanmean(t_sec_ORGNO3SO4NH4CHL_all,SO4_DL_all,t_now)
                NH4_write = cal_nanmean(t_sec_ORGNO3SO4NH4CHL_all,NH4_all,t_now)
                NH4_PREC_write = cal_nanmean(t_sec_ORGNO3SO4NH4CHL_all,NH4_PREC_all,t_now)
                NH4_DL_write = cal_nanmean(t_sec_ORGNO3SO4NH4CHL_all,NH4_DL_all,t_now)
                CHL_write = cal_nanmean(t_sec_ORGNO3SO4NH4CHL_all,CHL_all,t_now)
                CHL_PREC_write = cal_nanmean(t_sec_ORGNO3SO4NH4CHL_all,CHL_PREC_all,t_now)
                CHL_DL_write = cal_nanmean(t_sec_ORGNO3SO4NH4CHL_all,CHL_DL_all,t_now)
                # PAN
                LAT_PAN_write = cal_nanmean(t_sec_PAN_all,LAT_PAN_all,t_now) 
                LON_PAN_write = cal_nanmean(t_sec_PAN_all,LON_PAN_all,t_now)
                ALT_PAN_write = cal_nanmean(t_sec_PAN_all,ALT_PAN_all,t_now)
                THETA_PAN_write = cal_nanmean(t_sec_PAN_all,THETA_PAN_all,t_now)
                PAN_PPB_write = cal_nanmean(t_sec_PAN_all,PAN_PPB_all,t_now) 
                # VOC
                FOR_write = cal_nanmean(t_sec_VOC_all, FOR_all, t_now)
                FOR_PREC_write = cal_nanmean(t_sec_VOC_all, FOR_PREC_all, t_now)
                FOR_LOD_write = cal_nanmean(t_sec_VOC_all, FOR_LOD_all, t_now)
                MET_write = cal_nanmean(t_sec_VOC_all, MET_all, t_now)
                MET_PREC_write = cal_nanmean(t_sec_VOC_all, MET_PREC_all, t_now)
                MET_LOD_write = cal_nanmean(t_sec_VOC_all, MET_LOD_all, t_now)
                ACN_write = cal_nanmean(t_sec_VOC_all, ACN_all, t_now)
                ACN_PREC_write = cal_nanmean(t_sec_VOC_all, ACN_PREC_all, t_now)
                ACA_write = cal_nanmean(t_sec_VOC_all, ACA_all, t_now)
                ACA_PREC_write = cal_nanmean(t_sec_VOC_all, ACA_PREC_all, t_now)
                ACE_write = cal_nanmean(t_sec_VOC_all, ACE_all, t_now)
                ACE_PREC_write = cal_nanmean(t_sec_VOC_all, ACE_PREC_all, t_now)
                ISO_write = cal_nanmean(t_sec_VOC_all, ISO_all, t_now)
                ISO_PREC_write = cal_nanmean(t_sec_VOC_all, ISO_PREC_all, t_now)
                BEN_write = cal_nanmean(t_sec_VOC_all, BEN_all, t_now)
                BEN_PREC_write = cal_nanmean(t_sec_VOC_all, BEN_PREC_all, t_now)
                TOL_write = cal_nanmean(t_sec_VOC_all, TOL_all, t_now)
                TOL_PREC_write = cal_nanmean(t_sec_VOC_all, TOL_PREC_all, t_now)
                XYL_write = cal_nanmean(t_sec_VOC_all, XYL_all, t_now)
                XYL_PREC_write = cal_nanmean(t_sec_VOC_all, XYL_PREC_all, t_now)
                DMS_write = cal_nanmean(t_sec_VOC_all, DMS_all, t_now)
                DMS_PREC_write = cal_nanmean(t_sec_VOC_all, DMS_PREC_all, t_now)
                DMS_LOD_write = cal_nanmean(t_sec_VOC_all, DMS_LOD_all, t_now)
                MVK_write = cal_nanmean(t_sec_VOC_all, MVK_all, t_now)
                MVK_PREC_write = cal_nanmean(t_sec_VOC_all, MVK_PREC_all, t_now)
                MEK_write = cal_nanmean(t_sec_VOC_all, MEK_all, t_now)
                MEK_PREC_write = cal_nanmean(t_sec_VOC_all, MEK_PREC_all, t_now)
                # NO2_MINIDOAS
                NO2_MINIDOAS_write = cal_nanmean(t_sec_NO2_MINIDOAS_all, NO2_MINIDOAS_all, t_now)
                NO2_MINIDOAS_ERROR_write = cal_nanmean(t_sec_NO2_MINIDOAS_all, NO2_MINIDOAS_ERROR_all, t_now)

                write_text = [date_nc[j],strftime("%H:%M:%S", localtime(t_sec[j])), t_sec[j],\
                             lat[j],lon[j],H[j],H_round_to_hundred[j], u_ma[j], v_ma[j], w_ma[j],\
                             p_H[j], static_P[j], dynamic_P[j], TAT[j], TD[j], THETA_aircraft[j],\
                             THETA_V[j], TV[j], TS[j], RELHUM[j], MIXRATIOV_H2O[j], MIXRATIOM_H2O[j],\
                             SO2_write, CO_write, NO_write, NOy_write, O3_AMTEX_write, O3_FAIROCI_write,\
                             CH4_write, CO2_write, ORG_write, ORG_PREC_write, ORG_DL_write,\
                             NO3_write, NO3_PREC_write, NO3_DL_write, SO4_write, SO4_PREC_write, SO4_DL_write,\
                             NH4_write, NH4_PREC_write, NH4_DL_write, CHL_write, CHL_PREC_write, CHL_DL_write,\
                             LAT_PAN_write, LON_PAN_write, ALT_PAN_write, THETA_PAN_write, PAN_PPB_write,\
                             FOR_write, FOR_PREC_write, FOR_LOD_write, MET_write, MET_PREC_write, MET_LOD_write,\
                             ACN_write, ACN_PREC_write, ACA_write, ACA_PREC_write, ACE_write, ACE_PREC_write,\
                             ISO_write, ISO_PREC_write, BEN_write, BEN_PREC_write, TOL_write, TOL_PREC_write,\
                             XYL_write, XYL_PREC_write, DMS_write, DMS_PREC_write, DMS_LOD_write,\
                             MVK_write, MVK_PREC_write, MEK_write, MEK_PREC_write,\
                             NO2_MINIDOAS_write, NO2_MINIDOAS_ERROR_write]
                writer.writerow(write_text)
    fp.close()
    print('All files are finished!!!')


if __name__ == "__main__":
    time_start = time.time()
    main()
    print('All done!')
    print('Total time use: %.1f min' %((time.time() - time_start)/60))