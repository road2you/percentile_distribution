# 0) 코드 요약 ----------------------------------------------------------------
# 
# 백분위 기반 분포 계산 코드
# 출처 : 이준영, 박진서 (2019) 
#        학술논문 데이터로 바라 본 글로벌 과학기술연구 수준의 거시적 변동
#        (KISTI DATA INSIGHT 제7호)
#        http://mirian.kisti.re.kr/insight/insight.jsp
# 본 코드의 백분위 구간 비중 분할 방식은 아래 문헌의 알고리즘을 구현한 것임. 
# Waltman, L., & Schreiber, M. (2013). On the calculation of percentile-based 
# bibliometric indicators. Journal of the American Society for Information 
# Science and Technology, 64(2), 372–379. https://doi.org/10.1002/asi.22775
# 

# 데이터 처리에 활용한 패키지
library(data.table)
library(purrr)
library(stringr)

# 데이터 시각화에 활용한 패키지 
library(ggplot2)
library(ggthemes)
library(scales)
library(extrafont) 
library(wesanderson)
library(RColorBrewer)



# 1) 샘플데이터 반입  ------------------------------------------------------------

# OECD 기준 39개 분류 중 아래 데이터는 
# Nano-technology, Chemical Sciences 2개 분야만 수록

oecd_minor_ref_dt_sample <- fread("../data/oecd_minor_ref_dt_sample.txt")

oecd_minor_ref_dt_sample

# 데이터 구조 : 문헌고유번호 - 분류 - 년도 - 피인용수 
#
#                    uid        oecd_minor pubyear  tc
# 1: WOS:000206175400005 Chemical Sciences    1990   0
# 2: WOS:000207572600001 Chemical Sciences    1990   0
# 3: WOS:000207572600004 Chemical Sciences    1990   0
# 4: WOS:000207572600007 Chemical Sciences    1990   0
# 5: WOS:000207572600011 Chemical Sciences    1990   0



# 2) 백분위 계산용 기초 데이터 계산 ----------------------------------------------------

# 1) 분류코드 & pubyear 기준으로 citation 빈도 구하기 (Waltman's ci)

oecd_minor_frq_dt <-
  oecd_minor_ref_dt_sample[, .N,
                           by = .(oecd_minor, pubyear, tc)][order(oecd_minor, pubyear, tc)]

# 2) 분류코드 & pubyear 기준으로 citation 빈도의 합 구하기

oecd_minor_frq_dt[,
                  all_N := sum(N), by = .(oecd_minor, pubyear)]

# 3) 분류코드 & pubyear 기준으로 citation 누적 빈도 구하기

oecd_minor_frq_dt[,
                  cum_N := cumsum(N), by = .(oecd_minor, pubyear)]

# 4) 분류코드 & pubyear 기준으로 citation i 아래의 비율 (waltman qi)

oecd_minor_frq_dt[,
                  qi := unlist(map2(N, cum_N, function(x, y)
                    (y - x) / unique(all_N))), by = .(oecd_minor, pubyear)]


# 5) 분류코드 & pubyear 기준으로 citation i 아래의 비율 (waltman qi+1)

oecd_minor_frq_dt[, qi_plus := c(qi[-1], 1), by = .(oecd_minor, pubyear)]


# 6) 분류코드 & pubyear 기준으로 [qi, qi+1] 구간 구하기

oecd_minor_frq_dt[, q_range := (qi_plus - qi), by = .(oecd_minor, pubyear)]

oecd_minor_frq_dt 

# 최종 기초 데이터 계산 결과 형태
# 
#           oecd_minor pubyear  tc    N all_N cum_N        qi   qi_plus      q_range
# 1: Chemical Sciences    1990   0 7896 68883  7896 0.0000000 0.1146292 1.146292e-01
# 2: Chemical Sciences    1990   1 4811 68883 12707 0.1146292 0.1844722 6.984307e-02
# 3: Chemical Sciences    1990   2 4102 68883 16809 0.1844722 0.2440225 5.955025e-02
# 4: Chemical Sciences    1990   3 3533 68883 20342 0.2440225 0.2953123 5.128987e-02
# 5: Chemical Sciences    1990   4 3093 68883 23435 0.2953123 0.3402146 4.490223e-02
# 



# 3) 백분위 값과 분위 구간별 비중 계산 --------------------------------------------------

# 계산 병렬처리 효율을 위해 list 형태로 분할 

oecd_minor_dt_list <- split(oecd_minor_frq_dt, by = c("oecd_minor", "pubyear"))


# 각 분류별로 백분위 계산


# 백분위 구간 비율의 참조 테이블 생성 함수 
generate_percentile_ref <- function(k=100, pk, sk) {
   k_vec <- 1:k
   pk_vec <- c(0, pk)
   sk_vec <- sk
   pk_range_mat <- cbind(v1 = pk_vec[-length(pk_vec)], v2 = pk_vec[-1])
  return(list(kvec = k_vec, pk_mat = pk_range_mat, score_vec = sk_vec))
}


# 아래 percentile_score 함수의 인자로 투입할 참조 프레임 생성

# 100분위 구간 (default)
 R_100 <- generate_percentile_ref(k = 100, pk = seq(0.01, 1.00, 0.01), sk = 1:100)

# 10분위 구간
# R_10 <-  generate_percentile_ref(k = 10, pk = c(0.1, 0.2, 0.3, 0.4, 0.5, 
#                                                 0.6, 0.7, 0.8, 0.9, 1.00), 
#                        sk = 1:10)

# NSF Science & Engineering Indicator (6분위 구간)
# R_6 <- generate_percentile_ref(k = 6, pk = c(0.5, 0.75, 0.90, 0.95, 0.99, 1.00), 
#                        sk = 1:6)


# 백분위 점수 계산 함수 

percentile_score <-
  function(ori_cit, qi, qi_plus, q_range, perc_mat) {
    kvec <- perc_mat$kvec
    pk_mat <- perc_mat$pk_mat
    score_vec <- perc_mat$score_vec
    perc_score <- vector("numeric", length(qi))
    width_k <- floor(log(length(kvec), 10)) + 1
    perc_prop_var <-
      paste0("P", str_pad(1:length(kvec), width = width_k, pad = "0"))
    temp_score <- 0
    for (i in seq_along(perc_score)) {
      qi_temp <- qi[i]
      qi_plus_temp <- qi_plus[i]
      q_range_temp <- q_range[i]
      temp_dt <- data.table(ori_cit = ori_cit, perc_score = 0)
      temp_dt[, (perc_prop_var) := 0]
      for (j in seq_along(kvec)) {
        pk_minus1 <- pk_mat[j, 1]
        pk <- pk_mat[j, 2]
        s_temp <- score_vec[j]
        oik <- max(min(pk, qi_plus_temp) - max(pk_minus1, qi_temp), 0)
        set(
          temp_dt,
          i = i,
          j = 2L + j,
          value = oik / (q_range_temp)
        )
        temp_score <- oik / (q_range_temp) * s_temp + temp_score
      }
      set(temp_dt,
          i = i,
          j = 2L,
          value = temp_score)
    }
    return(temp_dt)
  }

# 최종 백분위 구간별 비율 계산 

oecd_percentile_res <-
  lapply(oecd_minor_dt_list, function(x)
    pmap_df(x[, .(tc, qi, qi_plus, q_range)],
            function(tc, qi, qi_plus, q_range)
              percentile_score(tc, qi, qi_plus, q_range, R_100)))

# 분야별 결과값 통합 
oecd_percentile_red_dt <- rbindlist(oecd_percentile_res)

oecd_percentile_red_dt[, oecd_minor := oecd_minor_frq_dt[, oecd_minor]]
oecd_percentile_red_dt[, pubyear := oecd_minor_frq_dt[, pubyear]]


oecd_merge_dt <- merge(oecd_minor_frq_dt, oecd_percentile_red_dt, 
                       by.x = c("oecd_minor", "pubyear", "tc"),
                       by.y = c("oecd_minor", "pubyear", "ori_cit"),
                       all.x = TRUE)


# 최종 계산 결과 형태 (P100까지의 칼럼 존재하나 일부만 남겼음)
# perc_score는 해당 피인용수(tc)가 걸쳐 있는 백분위 구간의 비중을 
# 고려하여 산출된 가중 평균 백분위 수임. 

#           oecd_minor pubyear  tc    N all_N cum_N        qi   qi_plus      q_range perc_score
# 1: Chemical Sciences    1990   0 7896 68883  7896 0.0000000 0.1146292 1.146292e-01   6.242302
# 2: Chemical Sciences    1990   1 4811 68883 12707 0.1146292 0.1844722 6.984307e-02  15.454968
# 3: Chemical Sciences    1990   2 4102 68883 16809 0.1844722 0.2440225 5.955025e-02  21.924166
# 4: Chemical Sciences    1990   3 3533 68883 20342 0.2440225 0.2953123 5.128987e-02  27.467577
# 5: Chemical Sciences    1990   4 3093 68883 23435 0.2953123 0.3402146 4.490223e-02  32.250954
#          P001       P002       P003       P004       P005       P006       P007       P008
# 1: 0.08723784 0.08723784 0.08723784 0.08723784 0.08723784 0.08723784 0.08723784 0.08723784
# 2: 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
# 3: 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
# 4: 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
# 5: 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000



# 만약 위에서 R_100 (default)으로 계산했다면, 
# 10분위 (이준영 & 박진서 (2019)) 방식으로 구간을 통합하여 분석하려면 
# 아래와 같은 묶음 작업 진행 

oecd_merge_dt2 <- oecd_merge_dt[, R1 := rowSums(.SD), .SDcols = 11:20]
oecd_merge_dt2 <- oecd_merge_dt2[, R2 := rowSums(.SD), .SDcols = 21:30]
oecd_merge_dt2 <- oecd_merge_dt2[, R3 := rowSums(.SD), .SDcols = 31:40]
oecd_merge_dt2 <- oecd_merge_dt2[, R4 := rowSums(.SD), .SDcols = 41:50]
oecd_merge_dt2 <- oecd_merge_dt2[, R5 := rowSums(.SD), .SDcols = 51:60]
oecd_merge_dt2 <- oecd_merge_dt2[, R6 := rowSums(.SD), .SDcols = 61:70]
oecd_merge_dt2 <- oecd_merge_dt2[, R7 := rowSums(.SD), .SDcols = 71:80]
oecd_merge_dt2 <- oecd_merge_dt2[, R8 := rowSums(.SD), .SDcols = 81:90]
oecd_merge_dt2 <- oecd_merge_dt2[, R9 := rowSums(.SD), .SDcols = 91:100]
oecd_merge_dt2 <- oecd_merge_dt2[, R10 := rowSums(.SD), .SDcols = 101:110]

oecd_R10_dt <- oecd_merge_dt2[, .(oecd_minor, pubyear, tc, perc_score,
                                  R1, R2, R3, R4, R5,
                                  R6, R7, R8, R9, R10)]
oecd_R10_dt

#           oecd_minor pubyear  tc perc_score        R1        R2        R3       R4 R5 R6 R7 R8 R9
# 1: Chemical Sciences    1990   0   6.242302 0.8723784 0.1276216 0.0000000 0.000000  0  0  0  0  0
# 2: Chemical Sciences    1990   1  15.454968 0.0000000 1.0000000 0.0000000 0.000000  0  0  0  0  0
# 3: Chemical Sciences    1990   2  21.924166 0.0000000 0.2607509 0.7392491 0.000000  0  0  0  0  0
# 4: Chemical Sciences    1990   3  27.467577 0.0000000 0.0000000 1.0000000 0.000000  0  0  0  0  0
# 5: Chemical Sciences    1990   4  32.250954 0.0000000 0.0000000 0.1043970 0.895603  0  0  0  0  0
#    R10
# 1:   0
# 2:   0
# 3:   0
# 4:   0
# 5:   0




# 4) 결과 분석 및 시각화 작업을 위한 데이터 작업 --------------------------------------------

# 위 결과데이터에 대상 문헌집합 매칭 
# 백분위수가 년도별 및 분야별 기준으로 설정한 영역에서 산출된 것이므로, 
# 피인용수 기준으로 매칭하며, 
# 문헌별 국가 정보를 )의 피인용수와 매칭
# 
# 이준영 & 박진서 (2019)에 수록된 8개 분야 중 
# 위에서 샘플 분야로 선택한 2개 분야, 
# Chemical Sciences, Nano-technology 로딩 

dt2minor <- fread("./data/dt2minor_sample.txt")


# 시각화를 위한 한글 글꼴 로딩  
font_import(pattern = 'NanumBarun') 
loadfonts(device = "win")
windowsFonts()



# 5) 최상위 10% 엑셀런스 및 최하위 10% 역엑셀런스 시계열 분석  ---------------------------------

# 데이터를 한국에 대해서만 추출
dt2minor_kr <- dt2minor[country2=="South Korea"]

# bottom_10, top_10의 문헌별 비중 값을 합산 

dt2minor_kr_both10 <- dt2minor_kr[,lapply(.SD[, c(9,18)], sum) , by = .(oecd_minor, pubyear)]

# 년도별 분석대상 문헌 수 계산하여 데이터 병합
dt2minor_kr_both10_ref <- dt2minor_kr[,.(n_article = length(unique(uid))),
                                      by = .(oecd_minor, pubyear)]
dt2minor_kr_both10_calc <- merge(dt2minor_kr_both10, dt2minor_kr_both10_ref,
                                 by.x = c("oecd_minor", "pubyear"),
                                 by.y = c("oecd_minor", "pubyear"),
                                 all.x = TRUE)

# bottm_10, top_10 10분위 구간에서 한국의 해당 분야 전체 논문에서 
# 각 구간에 속하는 문헌의 비율 (엑셀런스) 계산

dt2minor_kr_both10_res <- 
  dt2minor_kr_both10_calc[,
                          lapply(.SD[, c(1, 2)], function(x)
                            x / n_article * 100), by = .(oecd_minor, pubyear)]

# ggplot 패키지에 용이한 long 타입 형태의 데이터로 변환 
dt2minor_kr_both10_res_long <- melt(dt2minor_kr_both10_res, id.vars = 1:2, measure = 3:4,
                                    variable.name = "both10",
                                    value.name = "prop")

# 순서 정렬 
dt2minor_kr_both10_res_long <- 
  dt2minor_kr_both10_res_long[order(oecd_minor, both10, pubyear)]

dt2minor_kr_both10_res_long

# 최종 데이터 형태 
# 
#           oecd_minor pubyear both10      prop
# 1: Chemical Sciences    1990     R1  6.568496
# 2: Chemical Sciences    1991     R1  9.861503
# 3: Chemical Sciences    1992     R1  8.381822
# 4: Chemical Sciences    1993     R1 10.150435
# 5: Chemical Sciences    1994     R1  9.107963

kr_both10_chemical <- 
  ggplot(data = dt2minor_kr_both10_res_long[oecd_minor == "Chemical Sciences"], 
                    aes(x = pubyear, y = prop, color = both10)) +
  geom_hline(yintercept = 10, color = "royal blue") +
  theme_bw() +
  geom_line() +
  geom_smooth()  +
  ggtitle("한국의 인용 영향력 상/하위 비율 변화 (화학 분야)") +  
  xlab("출판년도") + ylab("비율") +
  scale_color_discrete(name = '10분위 구간 비율', labels = c("R1(하위10%)", "R10(상위10%)")) +
  scale_x_continuous(breaks = 1990:2018) + 
    theme(axis.text.x = element_text(angle = 60, hjust = 0.9, vjust = 0.9, size = 10, 
                                    family = "NanumBarunGothic"),
        axis.text.y = element_text(family = "NanumBarunGothic", size= 10),
        axis.title = element_text(vjust = 1, family = "NanumBarunGothic", size= 11),
        strip.text.x = element_text(size = 11, face = "bold", family = "NanumBarunGothic"),
        plot.title = element_text(size = 12, face = "bold", family = "NanumBarunGothic"),
        legend.title = element_text(size = 11, family = "NanumBarunGothic"),
        legend.text = element_text(size = 10, family = "NanumBarunGothic")
        )  

# png 형태로 그림 저장 
ggsave(plot = kr_both10_chemical, filename = "korea_both10_change_chemical.png", 
       # type = "cairo", 
       units = "in", width = 14, height = 8,    dpi = 1200)



# 6) 각 분야별 1990년대, 2000년대, 2010년대 분포 변화 -----------------------------------

# 대상 국가 
cntry_labs <-  c("영국", "독일", "일본", "중국", "대한민국", "미국")
names(cntry_labs) <- c("GBR", "Germany", "Japan", "Peoples R China", "South Korea", "USA")

# 10년 단위 분포로 변수 정의 
dt2minor[pubyear %in% 1990:1999, generation := "1990-1999"]
dt2minor[pubyear %in% 2000:2009, generation := "2000-2009"]
dt2minor[pubyear %in% 2010:2018, generation := "2010-2018"]

# 10년 단위 시대 + 해당 분야에서 각 국가별 문헌수 참조 테이블 생성 
dt2minor_generation_ref <- 
  unique(dt2minor[, .(n_by_minor_cntry2), by = .(country2, oecd_minor, generation)])
dt2minor_generation_ref <- 
  dt2minor_generation_ref[, lapply(.SD[,1], sum), by = .(country2, oecd_minor, generation)]

# 10년 단위 시대 + 해당 분야에서 각 국가별 10분위 구간의 값 합산 
dt2minor_generation_sum <- 
  dt2minor[, lapply(.SD[,9:18 ], sum), by = .(country2, oecd_minor, generation)]

# 분석대상 국가에 한정하여 데이터 추출 
dt2minor_generation_sum_cntry6 <- 
  dt2minor_generation_sum[country2 %in% names(cntry_labs)]

# ggplot 패키지 이용을 위한 long 타입 형태로 변환 
dt2minor_generation_sum_cntry6_long <- melt(dt2minor_generation_sum_cntry6,
                                            id.vars = 1:3, measure.vars = 4:13,
                                            variable.name = "P10", value.name = "prop")

dt2minor_generation_sum_cntry6_long

# 엑셀런스 계산을 위해 앞서 계산한 참조값 테이블과 병합
dt2minor_generation_merge <- merge(dt2minor_generation_sum_cntry6_long, dt2minor_generation_ref, 
                                   by.x = c("oecd_minor", "country2", "generation"), 
                                   by.y = c("oecd_minor", "country2", "generation"),
                                   all.x = TRUE)
# 엑셀런스 계산 
dt2minor_generation_merge[, prop2 := prop/n_by_minor_cntry2*100]

# 순서 정렬
dt2minor_generation_merge <- 
  dt2minor_generation_merge[order(oecd_minor, country2, generation, P10)]

# 최종 데이터 형태 
# 
#          oecd_minor country2 generation P10     prop n_by_minor_cntry2     prop2
# 1: Chemical Sciences      GBR  1990-1999  R1 3143.801             58909  5.336708
# 2: Chemical Sciences      GBR  1990-1999  R2 4060.087             58909  6.892134
# 3: Chemical Sciences      GBR  1990-1999  R3 4893.249             58909  8.306454
# 4: Chemical Sciences      GBR  1990-1999  R4 5592.752             58909  9.493883
# 5: Chemical Sciences      GBR  1990-1999  R5 6147.806             58909 10.436106


# Chemical Sciences 분야에 대해 분석

distribution_chemical <- ggplot(dt2minor_generation_merge[oecd_minor == "Chemical Sciences"], 
                   aes(x = P10, y = prop2, fill = country2)) +
  geom_bar(stat = "identity") + theme_bw() + 
  facet_grid(country2~generation,  labeller = labeller(country2 = cntry_labs)) + 
  geom_hline(yintercept = 10,
             color = "blue") +
  xlab("10분위 구간") + ylab("비율") +
  ggtitle("화학 분야에서 10분위 분포 변화") +
  scale_fill_discrete(name = '국가', 
                       labels = c("영국", "독일", "일본", "중국", "대한민국", "미국")) +
     theme(
       axis.text.x = element_text(angle = 60, hjust = 0.9, vjust = 0.9, size = 10, 
                                    family = "NanumBarunGothic"),
        axis.text.y = element_text(family = "NanumBarunGothic", size= 10),
        axis.title = element_text(vjust = 1, family = "NanumBarunGothic", size= 11),
        strip.text.x = element_text(size = 11, face = "bold", family = "NanumBarunGothic"),
        strip.text.y = element_text(size = 11, face = "bold", family = "NanumBarunGothic"),
        plot.title = element_text(size = 12, face = "bold", family = "NanumBarunGothic"),
        legend.title = element_text(size = 11, family = "NanumBarunGothic"),
        legend.text = element_text(size = 10, family = "NanumBarunGothic")) 

distribution_chemical

ggsave(plot = distribution_chemical, filename = "distribution_change_chemical.png", 
       # type = "cairo", 
       units = "in", width = 12, height = 14,    dpi = 1200)
