# Study 13 raw data

Place the SPSS data file here, renamed exactly:

```text
data/raw/study_13/DATA_Cleaned_and_coded_for_condition.sav
```

Required columns:

```text
condition
subj_total
obj_total
```

The uploaded SPSS syntax uses:

```spss
T-TEST GROUPS=condition(0 1)
  /VARIABLES=subj_total biasaware_total intmotiv_total behint_total belong_total obj_total
  /ES DISPLAY(TRUE).
```

The reproduction script treats `condition == 0` as control and `condition == 1` as education.
