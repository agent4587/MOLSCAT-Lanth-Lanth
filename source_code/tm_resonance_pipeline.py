#!/usr/bin/env python3
"""
tm_resonance_pipeline.py
=========================

Автоматизация поиска и характеризации резонансов Фешбаха для molscat-Tm2
(Этап 0 плана, см. CONTEXT_Tm_CPL9_v4.md, §2).

Запускается ЛОКАЛЬНО пользователем (нужен скомпилированный molscat-Tm2.exe
и рабочий каталог с pot-Tm2.f/base9-Tm_Tm_AI.f, уже собранные в бинарник).
Работает через subprocess, ничего не знает про физику сверх того, что
записано в шаблонах &INPUT/&BASIS9/&POTL.

------------------------------------------------------------------------
О КАЛИБРОВКЕ ПАРСЕРА
------------------------------------------------------------------------
Регэкспы парсера (`parse_coarse_scan`, `parse_ifconv_block`) взяты НЕ из
предположений, а напрямую из реальных исходников MOLSCAT 2025.0
(source_code/calca.f -- FORMAT 605/604 для таблицы сканирования;
source_code/locpol.f -- FORMAT 130/106/150/125 для характеризации
резонанса; source_code/mol.driver.f -- связь IFCONV=IDECAY+1), плюс
сверены с реальным образцом `molscat-Tm2_smoketest_v2.out` из репозитория
проекта. Тем не менее это НОВЫЙ код, ни разу не гонявшийся на настоящем
molscat-Tm2.exe (в песочнице, где он писался, нет ни бинарника, ни сети) --
перед систематическими прогонами (Этап 5-6) всё равно стоит:

    1. Прогнать coarse scan на маленьком диапазоне (пару полей).
    2. Вызвать: python tm_resonance_pipeline.py selftest-parse <файл.out> [ifconv]
    3. Если регэкспы не сработали (в реальности могут отличаться пробелы/
       переносы строк, которые сложно предсказать без живого вывода) --
       прислать мне фрагмент .out целиком, поправлю за один проход.
------------------------------------------------------------------------

Конвейер (соответствует 6 пунктам Этапа 0 в CONTEXT_Tm_CPL9_v4.md):
    1. generate_coarse_scan_input()   — input для грубого скана (LSCAN)
    2. parse_coarse_scan()            — .out -> DataFrame (B, Re(a), Im(a))
    3. detect_candidates()            — авто-детекция кандидатов
    4. generate_ifconv_input()        — input для точечной характеризации
    5. characterize_candidate()       — авто-эскалация IFCONV 1->2->3
    6. build_summary_table()          — финальный CSV с метаданными прогона
"""

from __future__ import annotations

import argparse
import time
import dataclasses
import math
import re
import subprocess
import sys
from pathlib import Path
from typing import Optional

import pandas as pd
import numpy as np


# =========================================================================
# 0. КОНФИГУРАЦИЯ / МЕТАДАННЫЕ ПРОГОНА
# =========================================================================

@dataclasses.dataclass
class RunMetadata:
    """Какие члены разложения V(R) = Σ_{k,i} V^(i)_k(R)·T^(i)_k реально
    включены в этот прогон — записывается в CSV для последующего сравнения
    статистики между разными комбинациями (Этап 5-6).

    ВАЖНО: используем ТОЛЬКО точные (k,i)-метки из Table 1, Tiesinga et al.,
    NJP 23, 085007 (2021) -- НЕ словесные названия ("dipole-dipole",
    "anisotropic VdW"), т.к. словесные названия оказались неоднозначны (см.
    переписку: было неясно, относится ли "anisotropic VdW" к уже
    реализованному (k=2,i=1), или к чему-то ещё не реализованному).
    Таблица 1 для Tm2 (k, i, оператор):
      (0,1)=I изотропный;
      (2,1)=[j1⊗j1]2+[j2⊗j2]2 (ранг-2, Таблица 5);
      (0,2)=[j1⊗j2]0;
      (2,2)=[j1⊗j2]2 (ЭТО магнитный диполь-диполь на больших R:
             V^(2)_2(R) -> D^(2)_2/R^3 + C^(2)_2/R^6);
      (0,3)=[[j1⊗j1]2⊗[j2⊗j2]2]0;
      (2,3)=[[j1⊗j1]2⊗[j2⊗j2]2]2;
      (4,1)=[[j1⊗j1]2⊗[j2⊗j2]2]4 (ранг-4).
    В base9-Tm_Tm_AI.f реализованы и верифицированы (0,1) и (2,1) (VSTAR,
    pot-Tm2.f), а также D^(2)_2/R^3 часть (2,2) (аналитический член в
    &POTL, см. DIPBLK/TENSX в base9-Tm_Tm_AI.f) -- ЭТО магнитный
    диполь-диполь. Слабая дисперсионная добавка C^(2)_2/R^6 к тому же
    V^(2)_2(R) НЕ включена (на 2+ порядка меньше D^(2)_2/R^3 во всём
    диапазоне интереса)."""
    tensor_terms: str = "k0i1,k2i1,k2i2"  # список включённых (k,i)-членов
                                       # V(R), запятая между ними;
                                       # k2i2 подразумевает dipdip_scale=1
                                       # (см. Config.dipdip_scale) -- это
                                       # ТОЛЬКО метаданные для CSV, саму
                                       # физику переключает dipdip_scale
    hyperfine_zeeman: bool = True     # одноатомные члены (HFSPLA/GSA в
                                       # BAS9IN) -- НЕ часть разложения V(R),
                                       # отдельная ось; сейчас всегда включены
                                       # (переключатель "выкл" -- Этап 4, ещё
                                       # не реализован)
    lmax: int = 4                     # LMAX БАЗИСА (число партиальных волн
                                       # в разложении по l), это НЕ тензорный
                                       # ранг k -- не путать при чтении CSV!
    dipdip_scale: float = 1.0         # ДУБЛИРУЕТ Config.dipdip_scale ТОЛЬКО
                                       # для записи в CSV (см. run_full_pipeline)
    v2_scale: float = 1.0             # ДУБЛИРУЕТ Config.v2_scale ТОЛЬКО
                                       # для записи в CSV (см. run_full_pipeline)
    label: str = "k0i1_k2i1_k2i2"


@dataclasses.dataclass
class Config:
    molscat_exe: Path
    work_dir: Path                   # cwd для запуска molscat-Tm2.exe (там,
                                       # где лежит бинарник и его зависимости)
    output_dir: Path                 # куда писать .input/.out (ОТДЕЛЬНО от
                                       # work_dir, чтобы не засорять source_code)
    coarse_template: Path            # шаблон &INPUT/&BASIS/&BASIS9/&POTL
    energy_K: float = 1.0e-6         # энергия столкновения (Кельвин), ~1 мкК
    jtot: int = -12
    ibfix: int = 2
    jstep: int = 2
    lmax: int = 4
    ured: float = 84.467109
    lambda_terms: str = "0, 1"       # MXLAM/LAMBDA как в molscat-Tm2_fieldscan_v2
    dipdip_scale: float = 1.0        # множитель перед D^(2)_2 (магнитный
                                       # диполь-диполь, k=2,i=2, [j1⊗j2]_2q --
                                       # см. DIPBLK/TENSX в base9-Tm_Tm_AI.f).
                                       # 1.0 = реальный физический коэффициент
                                       # (D2_DIPDIP_CM_ANG3 ниже), 0.0 = член
                                       # выключен (V^(2)_2(R)=0, как до его
                                       # добавления), любое другое значение --
                                       # искусственное масштабирование для
                                       # сравнительных прогонов.
    v2_scale: float = 1.0            # множитель перед всем коэффициентом
                                       # связи блока k=2,i=1 (анизотропный
                                       # V^(1)_2(R), Table 5 -- см. V2SCALE
                                       # в &BASIS9 / base9-Tm_Tm_AI.f). Физически
                                       # эквивалентно масштабированию V^(1)_2(R)
                                       # целиком для всех R (коэффициент связи и
                                       # радиальная функция входят произведением).
                                       # 1.0 = реальная физика (дефолт), 0.0 =
                                       # анизотропия k=2,i=1 выключена.
    ichan_guess: int = 1             # входной канал в базисе (см. "грабли" v4 §4)
    monqn: Optional[str] = None      # "2F1,2mF1,2F2,2mF2" (доубленные, через
                                       # запятую) -- если задано, EREF считается
                                       # через THRSH9 (Брейт-Рабин) по этим
                                       # квантовым числам вместо IREF=1 (см.
                                       # THRSH9 в base9-Tm_Tm_AI.f). ICHAN всё
                                       # равно нужен отдельно (см. ichan_guess/
                                       # candidate["channel"]) -- MONQN решает
                                       # только вырождение/выбор EREF, не выбор
                                       # канала распространения.


# =========================================================================
# 1. ГЕНЕРАЦИЯ INPUT ДЛЯ ГРУБОГО СКАНА
# =========================================================================

# D^(2)_2 = -sqrt(6)*alpha^2*(g_j/2)^2*E_h*a0^3 (Tiesinga et al., NJP 23,
# 085007 (2021), Sec 2.3), g_j=1.14119 (=GSA в base9-Tm_Tm_AI.f), переведено
# a0^3 -> Angstrom^3 (R в &POTL -- ангстремы, как и в pot-Tm2.f), теми же
# константами, что реально линкует Makefile (physical_constants_module.f,
# 2022 CODATA/BIPM: 1/alpha=137.035999177, E_h=2.1947463136314e5 cm^-1,
# bohr_to_Angstrom=0.529177210544). Умножается на Config.dipdip_scale и
# подставляется как единственный аналитический член (NPOWER=-3) 3-го блока
# &POTL -- см. DIPBLK/TENSX и POTIN9 в base9-Tm_Tm_AI.f.
D2_DIPDIP_CM_ANG3 = -1.38117889454613  # cm^-1 * Angstrom^3, при dipdip_scale=1

COARSE_TEMPLATE = """\
 &INPUT
    LABEL  =  '{label}',
    URED   =  {ured},
    IPRINT =   6,
    RMIN   =   3.0,  RMID   =  21.0,   RMAX   = 15.0E3, IRMSET = 0,
    IPROPS =   6,    DR     =   0.002,
    IPROPL =   9,    TOLHIL =   1.E-7,
    JTOTL  = {jtot},    JTOTU  = {jtot},     IBFIX  =  {ibfix},     JSTEP  = {jstep},
    EUNITS =   2,    NNRG   =   1,     ENERGY =  {energy},
                     DTOL   =   1.E-6,
    FLDMIN =    {fmin}, FLDMAX = {fmax},  DFIELD =   {dfield},
    IREF   =   {iref},    ICHAN  =   {ichan},{monqn_clause}
    LASTIN =   1,
 /

 &BASIS
    ITYPE  = 9,
 /

 &BASIS9
    LMAX   =  {lmax},   V2SCALE = {v2_scale},
 /

 &POTL
    MXLAM  = 3, LAMBDA =  {lambda_terms}, 2,
                NTERM  = -1, -1,  1,
                NPOWER =              -3,
                A      =              {dipdip_a},
 /
"""

IFCONV_TEMPLATE = """\
 &INPUT
    LABEL  =  '{label}',
    URED   =  {ured},
    IPRINT =   7,
    RMIN   =   3.0,  RMID   =  21.0,   RMAX   = 15.0E3, IRMSET = 0,
    IPROPS =   6,    DR     =   0.002,
    IPROPL =   9,    TOLHIL =   1.E-7,
    JTOTL  = {jtot},    JTOTU  = {jtot},     IBFIX  =  {ibfix},     JSTEP  = {jstep},
    EUNITS =   2,    NNRG   =   1,     ENERGY =  {energy},
                     DTOL   =   0.01,
    FLDMIN =  {fmin}, FLDMAX = {fmax},
    IREF   =   {iref},    ICHAN  =   {ichan},{monqn_clause}     IFCONV =   {ifconv}{iphsum_clause},
    LASTIN =   1,
 /

 &BASIS
    ITYPE  = 9,
 /

 &BASIS9
    LMAX   =  {lmax},   V2SCALE = {v2_scale},
 /

 &POTL
    MXLAM  = 3, LAMBDA =  {lambda_terms}, 2,
                NTERM  = -1, -1,  1,
                NPOWER =              -3,
                A      =              {dipdip_a},
 /
"""

# =========================================================================
# 1b. АЛЬТЕРНАТИВНЫЙ МЕТОД: ПОИСК РЕЗОНАНСОВ ЧЕРЕЗ FIELD (не molscat)
# =========================================================================
#
# field.exe (fld.driver.f) решает ОБРАТНУЮ задачу к bound: "при заданной
# целевой энергии найди ПОЛЕ, при котором связанный уровень её достигает"
# (см. fld.driver.f:21-22 -- "TO FIND THE APPLIED FIELD WHERE A BOUND
# STATE HAS A SPECIFIED ENERGY. WHEN ENERGY IS SET TO A THRESHOLD, THIS
# FINDS RESONANCE POSITIONS"). С ENERGY=0.0 (относительно порога, через
# тот же механизм MONQN/THRSH9, что и в bound) это даёт B0 НАПРЯМУЮ,
# бисекцией по полю -- без грубого скана molscat, без эвристик
# detect_candidates по скачку Re(a)/Im(a), и может найти НЕСКОЛЬКО
# резонансов в одном окне FLDMIN-FLDMAX за один прогон (программа сама
# считает разницу числа узлов на FLDMIN/FLDMAX и ищет ровно столько
# уровней, сколько эта разница подразумевает). Портировано из
# toy_resonance_pipeline.py (find-resonances) -- см. там же обсуждение
# про то, что FLDMIN/FLDMAX задают только СТАРТОВОЕ окно для field.exe/
# LOCPOL, а не жёсткую границу поиска.
#
# ВАЖНО: MONQN здесь ОБЯЗАТЕЛЕН (не опционален, как в COARSE/IFCONV
# TEMPLATE) -- без него не от чего считать "порог", и вся идея метода
# теряет смысл. IREF всегда 0 (жёстко, не namelist-параметр здесь). Для
# Tm2 (I=1/2) MONQN -- 4 доубленных числа '2F1,2mF1,2F2,2mF2', напр.
# '8,-8,8,-8' для полностью растянутого состояния F1=F2=4,mF1=mF2=-4 (см.
# ABS(2F-ISA)==1 в THRSH9 base9-Tm_Tm_AI.f -- это ДРУГАЯ семантика, чем в
# toy_resonance_pipeline.py, где I=0 и MONQN(1)=ISA ТОЧНО).
#
# RMATCH используется вместо RMID (терминология bound-state программ
# bound/field, а не молскатовского RMID) -- см. field-Rb2.input,
# fld.driver.f: "RMATCH IS THE MATCHING POINT" между исходящим и
# входящим решениями, тогда как RMID -- точка смены МЕТОДА
# распространения (для field/bound это не то же самое, что для
# рассеяния, но численно можно использовать те же значения, что и в
# COARSE_TEMPLATE -- проверено вживую в toy_resonance_pipeline.py).
FIELD_TEMPLATE = """\
 &INPUT
    LABEL  =  '{label}',
    URED   =  {ured},
    IPRINT =   6,
    RMIN   =   3.0,  RMATCH =  21.0,   RMAX   = 15.0E3, IRMSET = 0,
    IPROPS =   6,    DR     =   0.002,
    IPROPL =   9,    TOLHIL =   1.E-7,
    JTOTL  = {jtot},    JTOTU  = {jtot},     IBFIX  =  {ibfix},     JSTEP  = {jstep},
    EUNITS =   1,    NNRG   =   1,     ENERGY =  0.0,
                     DTOL   =   1.E-6,
    FLDMIN =  {fldmin}, FLDMAX = {fldmax},
    MONQN  =   {monqn_values},
    LASTIN =   1,
 /

 &BASIS
    ITYPE  = 9,
 /

 &BASIS9
    LMAX   =  {lmax},   V2SCALE = {v2_scale},
 /

 &POTL
    MXLAM  = 3, LAMBDA =  0, 1, 2,
                NTERM  = -1, -1,  1,
                NPOWER =              -3,
                A      =              {dipdip_a},
 /
"""


def _iref_and_monqn_clause(monqn: Optional[str]) -> tuple[int, str]:
    """Переводит cfg.monqn ("8,-8,8,-8" и т.п.) в пару (IREF, MONQN-клоза
    для namelist). Если monqn не задан -- старое поведение (IREF=1, без
    MONQN). Если задан -- IREF=0 (обязательное требование THRSH9, см.
    'ERROR - THRSH9 CALLED WITH POSITIVE IREF' в base9-Tm_Tm_AI.f) и
    добавляется строка MONQN = <значения>. Значения ДОЛЖНЫ быть уже
    доубленными (2*F, 2*mF), как того ждёт THRSH9 -- функция сама их не
    удваивает и не проверяет физическую осмысленность (это входные
    квантовые числа пользователя, см. THRSH9 для правил ABS(2F-ISA)==1
    и чётности 2F+2mF)."""
    if not monqn:
        return 1, ""
    values = [v.strip() for v in monqn.split(",")]
    if len(values) != 4:
        raise ValueError(
            f"--monqn должен содержать РОВНО 4 доубленных числа "
            f"(2*F1,2*mF1,2*F2,2*mF2), получено {len(values)}: {monqn!r}"
        )
    return 0, "\n    MONQN  =   " + ", ".join(values) + ","


def _fortran_e(x: float) -> str:
    """Форматирует число в стиле Fortran-литерала (1.E-6 вместо 1e-06) --
    численно эквивалентно, но точь-в-точь как в рабочих .input пользователя."""
    s = f"{x:.6E}"
    mantissa, exp = s.split("E")
    mantissa = mantissa.rstrip("0").rstrip(".") or "0"
    return f"{mantissa}E{exp}"


def _dipdip_a(cfg: Config) -> str:
    """Коэффициент A для 3-го блока &POTL (D^(2)_2/R^3, k=2,i=2), масштабированный
    Config.dipdip_scale. Фиксированный формат с плавающей точкой (не %e) --
    Fortran namelist READ одинаково понимает оба, но так удобнее читать
    сгенерированный .input глазами."""
    return f"{cfg.dipdip_scale * D2_DIPDIP_CM_ANG3:.14f}"


def generate_coarse_scan_input(cfg: Config, fmin: float, fmax: float,
                                 dfield: float, label: str,
                                 out_path: Path) -> Path:
    iref, monqn_clause = _iref_and_monqn_clause(cfg.monqn)
    text = COARSE_TEMPLATE.format(
        label=label, ured=cfg.ured, jtot=cfg.jtot, ibfix=cfg.ibfix,
        jstep=cfg.jstep, energy=_fortran_e(cfg.energy_K), fmin=fmin, fmax=fmax,
        dfield=dfield, ichan=cfg.ichan_guess, lmax=cfg.lmax,
        lambda_terms=cfg.lambda_terms, iref=iref, monqn_clause=monqn_clause,
        dipdip_a=_dipdip_a(cfg), v2_scale=cfg.v2_scale,
    )
    out_path.write_text(text)
    return out_path


def generate_ifconv_input(cfg: Config, fmin: float, fmax: float,
                            ichan: int, ifconv: int,
                            label: str, out_path: Path) -> Path:
    """Узкое бракетирующее окно вокруг кандидата, IFCONV=1/2/3/4.
    fmin/fmax — границы окна, внутри которого LOCPOL сам подбирает 3 точки
    сходимости (DFIELD тут НЕ участвует — в отличие от грубого скана,
    подтверждено реальным рабочим .input пользователя).

    ВАЖНО (найдено на реальном отказавшем прогоне, 260 Гс): для IFCONV=4
    режим характеризации (LCHAR) в mol.driver.f включается ТОЛЬКО при
    одновременном выполнении IFCONV=4 И IPHSUM>0 -- без IPHSUM MOLSCAT
    молча откатывается к обычному сканированию по сетке (никакого LOCPOL,
    никакой сходимости, просто печатает K-DEPENDENT... на нескольких полях
    и завершается без единого узнаваемого маркера, отсюда status='unknown').
    Юнит 9 для IPHSUM пишет доп. файл (fort.9) в cfg.work_dir -- он
    перезаписывается на каждый прогон, не хранит ничего критичного."""
    iphsum_clause = ",     IPHSUM =   9" if ifconv == 4 else ""
    iref, monqn_clause = _iref_and_monqn_clause(cfg.monqn)
    text = IFCONV_TEMPLATE.format(
        label=label, ured=cfg.ured, jtot=cfg.jtot, ibfix=cfg.ibfix,
        jstep=cfg.jstep, energy=_fortran_e(cfg.energy_K), fmin=fmin, fmax=fmax,
        ichan=ichan, ifconv=ifconv, lmax=cfg.lmax,
        lambda_terms=cfg.lambda_terms, iphsum_clause=iphsum_clause,
        iref=iref, monqn_clause=monqn_clause, dipdip_a=_dipdip_a(cfg),
        v2_scale=cfg.v2_scale,
    )
    out_path.write_text(text)
    return out_path


def generate_field_scan_input(cfg: Config, fldmin: float, fldmax: float,
                                monqn: str, label: str, out_path: Path) -> Path:
    """Input для field.exe (FIELD_TEMPLATE) -- ищет поле(я), где связанный
    уровень пересекает порог (ENERGY=0 относительно порога по MONQN), внутри
    [fldmin, fldmax]. См. комментарий перед FIELD_TEMPLATE.

    monqn ОБЯЗАТЕЛЕН здесь (в отличие от cfg.monqn для molscat-шаблонов) --
    без него не определён порог, относительно которого ищется пересечение."""
    values = [v.strip() for v in monqn.split(",")]
    if len(values) != 4:
        raise ValueError(
            f"monqn для find-resonances должен содержать РОВНО 4 доубленных "
            f"числа (2*F1,2*mF1,2*F2,2*mF2), получено {len(values)}: {monqn!r}"
        )
    text = FIELD_TEMPLATE.format(
        label=label, ured=cfg.ured, jtot=cfg.jtot, ibfix=cfg.ibfix,
        jstep=cfg.jstep, fldmin=fldmin, fldmax=fldmax,
        monqn_values=", ".join(values), lmax=cfg.lmax, v2_scale=cfg.v2_scale,
        dipdip_a=_dipdip_a(cfg),
    )
    out_path.write_text(text)
    return out_path


# =========================================================================
# 1c. ПОДСЧЁТ ЧИСЛА РЕЗОНАНСОВ В ОКНЕ (без поиска точных позиций)
# =========================================================================
#
# field.exe САМ считает число узлов волновой функции на FLDMIN (NODLO) и на
# FLDMAX (NODHI) -- это первое, что он делает, ДО дорогого поиска точных
# позиций (см. fld.driver.f:743-848, "DO IFLD=1,NFIELD" -- IFLD=1 -> FLD=
# FLDMIN, IFLD=2 -> FLD=FLDMAX, ровно 2 вычисления). Разница |NODHI-NODLO|
# -- это и есть число уровней, пересекающих порог между FLDMIN и FLDMAX
# (= число резонансов в окне), и печатается САМИМ field.exe как
# "NODE COUNT IN/DECREASES ... SEEK N STATE(S) IN THIS INTERVAL" (или
# "NO NODES BETWEEN FLDMAX AND FLDMIN", если резонансов нет вообще) --
# см. fld.driver.f:855-881.
#
# ВАЖНО: MXCALC (namelist &INPUT, дефолт 1000) -- максимальное число
# "тяжёлых" вычислений (распространений уравнений) за прогон; каждое
# вычисление NODLO/NODHI -- одно такое. Проверка `NCALC.GE.MXCALC` стоит
# ВНУТРИ цикла IFLD=1,2 (fld.driver.f:845) -- значит MXCALC=2 оборвал бы
# прогон СРАЗУ ПОСЛЕ второго вычисления (NODHI), но ДО печати нужного нам
# сообщения (которое печатается ПОСЛЕ конца цикла, fld.driver.f:861-881).
# MXCALC=3 даёт циклу нормально завершиться (оба NODLO/NODHI посчитаны,
# сообщение напечатано), и обрывает прогон максимум на ОДНОМ шаге
# дорогого поиска точных позиций -- вместо полного перебора всех N
# состояний, как в find-resonances. Проверено эмпирически (см. переписку
# про lmax=6, 0.1-20 Гс): экономит время именно на широких/густых окнах,
# где количество резонансов велико.
NODE_COUNT_TEMPLATE = """\
 &INPUT
    LABEL  =  '{label}',
    URED   =  {ured},
    IPRINT =   6,
    RMIN   =   3.0,  RMATCH =  21.0,   RMAX   = 15.0E3, IRMSET = 0,
    IPROPS =   6,    DR     =   0.002,
    IPROPL =   9,    TOLHIL =   1.E-7,
    JTOTL  = {jtot},    JTOTU  = {jtot},     IBFIX  =  {ibfix},     JSTEP  = {jstep},
    EUNITS =   1,    NNRG   =   1,     ENERGY =  0.0,
                     DTOL   =   1.E-6,
    FLDMIN =  {fldmin}, FLDMAX = {fldmax},
    MONQN  =   {monqn_values},
    MXCALC =   3,
    LASTIN =   1,
 /

 &BASIS
    ITYPE  = 9,
 /

 &BASIS9
    LMAX   =  {lmax},   V2SCALE = {v2_scale},
 /

 &POTL
    MXLAM  = 3, LAMBDA =  0, 1, 2,
                NTERM  = -1, -1,  1,
                NPOWER =              -3,
                A      =              {dipdip_a},
 /
"""


def generate_node_count_input(cfg: Config, fldmin: float, fldmax: float,
                                monqn: str, label: str, out_path: Path) -> Path:
    """Input для field.exe с MXCALC=3 -- см. комментарий перед
    NODE_COUNT_TEMPLATE. monqn ОБЯЗАТЕЛЕН, как и в generate_field_scan_input."""
    values = [v.strip() for v in monqn.split(",")]
    if len(values) != 4:
        raise ValueError(
            f"monqn для count-resonances должен содержать РОВНО 4 доубленных "
            f"числа (2*F1,2*mF1,2*F2,2*mF2), получено {len(values)}: {monqn!r}"
        )
    text = NODE_COUNT_TEMPLATE.format(
        label=label, ured=cfg.ured, jtot=cfg.jtot, ibfix=cfg.ibfix,
        jstep=cfg.jstep, fldmin=fldmin, fldmax=fldmax,
        monqn_values=", ".join(values), lmax=cfg.lmax, v2_scale=cfg.v2_scale,
        dipdip_a=_dipdip_a(cfg),
    )
    out_path.write_text(text)
    return out_path


# =========================================================================
# 2. ЗАПУСК MOLSCAT
# =========================================================================

def run_molscat(cfg: Config, input_path: Path) -> tuple[Path, float]:
    """Запускает molscat-Tm2.exe < input_path > input_path.with_suffix('.out').
    input_path/output .out живут в cfg.output_dir; сам процесс запускается с
    cwd=cfg.work_dir (там, где бинарник и его зависимости -- НЕ то же самое,
    что output_dir, см. Config).
    Возвращает (путь к .out, время выполнения в секундах). Не бросает
    исключение при ненулевом коде возврата — molscat иногда возвращает его
    даже при штатных сообщениях типа NOPEN CHANGED; решение о фатальности
    принимает вызывающий код по содержимому .out.

    БЕЗ таймаута -- ждёт завершения процесса сколько потребуется (см.
    переписку про count-resonances: жёсткий timeout=600 обрывал честно
    считающий field.exe на большом LMAX; вместо подбора "правильного"
    числа таймаут убран совсем)."""
    out_path = input_path.with_suffix(".out")
    t0 = time.perf_counter()
    with open(input_path, "r") as fin, open(out_path, "w") as fout:
        subprocess.run(
            [str(cfg.molscat_exe)],
            stdin=fin, stdout=fout, stderr=subprocess.STDOUT,
            cwd=cfg.work_dir,
        )
    elapsed = time.perf_counter() - t0
    return out_path, elapsed


def _fmt_duration(seconds: float) -> str:
    """Человекочитаемая длительность: 45s / 3m12s / 1h05m."""
    seconds = int(round(seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h{m:02d}m"
    if m:
        return f"{m}m{s:02d}s"
    return f"{s}s"


# =========================================================================
# 3. ПАРСИНГ ГРУБОГО СКАНА
# =========================================================================
#
# Формат подтверждён по реальным исходникам MOLSCAT (2025.0) в проекте:
#   source_code/calca.f (FORMAT 605/604) и реальному образцу
#   molscat-Tm2_smoketest_v2.out. Пример строки поля:
#     EFV SET     1:  MAGNETIC Z FIELD =   1000.000000     GAUSS
#   Пример заголовка таблицы:
#     K-DEPENDENT SCATTERING LENGTHS/VOLUMES/HYPERVOLUMES FOR CHANNELS...
#     CHAN   L POW     WVEC*ANGSTROM         RE(A)/ANGSTROM        IM(A)/ANGSTROM
#   Пример строки данных (FORMAT 604: 2I5,I4,4(ES21.13E3,2X), пишутся 3 из 4):
#     1    0   1 1.8661598530644E-003   1.8565979927580E+002  -4.1942512160504E-004
#   Столбцы: CHAN, L, POW, WVEC, RE(A), IM(A).

RE_EFV_SET = re.compile(
    r"EFV\s+SET\s+(\d+):.*?FIELD\s*=\s*([\-\d.]+(?:[eE][\-+]?\d+)?)",
    re.IGNORECASE,
)

RE_SCATLEN_HEADER = re.compile(r"K-DEPENDENT SCATTERING LENGTH", re.IGNORECASE)

_NUM = r"[\-\d.]+(?:[eE][\-+]?\d+)?"
RE_SCATLEN_ROW = re.compile(
    rf"^\s*(\d+)\s+(\d+)\s+(\d+)\s+({_NUM})\s+({_NUM})\s+({_NUM})\s*$"
)


def parse_coarse_scan(out_path: Path, target_ichan: Optional[int] = None) -> pd.DataFrame:
    """Возвращает DataFrame [field_G, channel, L, pow, wvec, re_a, im_a].
    Если target_ichan задан — фильтрует только этот канал (входной канал
    базиса, см. "грабли": ICHAN фиксирован в базисе, не "плывущая" строка
    таблицы открытых каналов; ВАЖНО: таблица печатает только НИЗКОЭНЕРГЕТИЧНЫЕ
    открытые каналы (DD<AWVMAX в calca.f) — target_ichan здесь фильтрует
    именно по этому CHAN-индексу, а не по позиции строки)."""
    text = out_path.read_text(errors="replace")
    lines = text.splitlines()

    rows = []
    current_field = None
    in_table = False

    for line in lines:
        m_field = RE_EFV_SET.search(line)
        if m_field:
            current_field = float(m_field.group(2))
            in_table = False
            continue

        if RE_SCATLEN_HEADER.search(line):
            in_table = True
            continue

        if in_table:
            m_row = RE_SCATLEN_ROW.match(line)
            if m_row:
                chan, L, pow_, wvec, re_a, im_a = m_row.groups()
                rows.append({
                    "field_G": current_field,
                    "channel": int(chan),
                    "L": int(L),
                    "pow": int(pow_),
                    "wvec": float(wvec),
                    "re_a": float(re_a),
                    "im_a": float(im_a),
                })
            elif line.strip() == "":
                in_table = False  # конец таблицы для этого EFV SET

    df = pd.DataFrame(rows)
    if df.empty:
        print(f"[WARN] Ничего не распарсено из {out_path}. "
              f"Похоже, IPRINT<6 (нет автопечати K-DEPENDENT...) либо "
              f"energía столкновения слишком высокая (DD>=AWVMAX — строка "
              f"вообще не печатается, см. calca.f). Запусти `selftest-parse` "
              f"и пришли фрагмент .out.",
              file=sys.stderr)
        return df

    if target_ichan is not None:
        df = df[df.channel == target_ichan].reset_index(drop=True)
    return df.sort_values("field_G").reset_index(drop=True)


# =========================================================================
# 3b. ПАРСИНГ FIELD (поиск резонансов напрямую, без грубого скана)
# =========================================================================
#
# Формат взят из ИСХОДНИКА подпрограммы, которая печатает эти строки --
# source_code/prnsum.f, SUBROUTINE PRCONV, ЧЕТЫРЕ FORMAT-а (100/200/300/400):
#   100: "CONVERGED ON STATE NUMBER I5 AT A = value unit"            (чисто)
#   200: та же строка + следующая "  BUT LAST NODE COUNT = I5 NOT AS DESIRED"
#        (нашли уровень, но число узлов не совпало с ожидаемым -- возможно,
#        случайное вырождение/пересечение уровней рядом)
#   300: та же строка + " BUT VARIABLE IS OUTSIDE RANGE" (сошлось за
#        пределами [FLDMIN,FLDMAX] -- обычно означает, что реального
#        пересечения в окне НЕТ, есть только около границы)
#   400: "NOT CONVERGED ON STATE NUMBER I5....) -- не сошлось за NITER
#        итераций, значение -- текущая грубая оценка, НЕ доверять как B0.
# Портировано из toy_resonance_pipeline.py, где проверено вживую на
# field-Toy.exe (2 резонанса в одном окне) -- формат общий для field.exe
# независимо от base9-модуля (это код prnsum.f, не base9-специфичный).
_RE_CONVERGED = re.compile(
    r"CONVERGED ON STATE NUMBER\s+(\d+)\s+AT\s+.+?=\s*"
    r"([-+]?\d+\.?\d*(?:[DdEe][-+]?\d+)?)\s+(\S+)"
)
_RE_NOT_CONVERGED = re.compile(r"NOT CONVERGED ON STATE NUMBER\s+(\d+)")
_RE_NO_NODES = re.compile(r"NO NODES BETWEEN FLDMAX AND FLDMIN")


def parse_field_resonances(out_path: Path) -> pd.DataFrame:
    """Извлекает найденные field.exe позиции резонансов (поля, где связанный
    уровень пересекает ENERGY=0 относительно порога) из .out файла.
    Возвращает DataFrame с колонками state, field_G, status, unit -- пустой
    (с этими же колонками), если резонансов в окне не было вообще (найдена
    строка 'NO NODES BETWEEN FLDMAX AND FLDMIN' -- ЭТО НЕ ошибка парсинга,
    а явный сигнал MOLSCAT "здесь ничего нет", в отличие от detect_candidates,
    где пустой результат неотличим от `IPRINT<6`/иных проблем без ручной
    проверки .out)."""
    lines = out_path.read_text(errors="replace").splitlines()
    rows: list[dict] = []
    found_no_nodes = False
    for i, line in enumerate(lines):
        if _RE_NO_NODES.search(line):
            found_no_nodes = True
            continue
        m_not = _RE_NOT_CONVERGED.search(line)
        if m_not:
            rows.append({"state": int(m_not.group(1)), "field_G": None,
                         "status": "not_converged", "unit": None})
            continue
        m = _RE_CONVERGED.search(line)
        if not m:
            continue
        state, value, unit = m.groups()
        status = "ok"
        if "OUTSIDE RANGE" in line:
            status = "out_of_range"
        elif i + 1 < len(lines) and "BUT LAST NODE COUNT" in lines[i + 1]:
            status = "node_mismatch"
        rows.append({
            "state": int(state),
            "field_G": float(value.replace("D", "E").replace("d", "e")),
            "status": status, "unit": unit,
        })
    df = pd.DataFrame(rows, columns=["state", "field_G", "status", "unit"])
    if df.empty and not found_no_nodes:
        print(f"[WARN] Ни одной строки CONVERGED/NOT CONVERGED/NO NODES не "
              f"найдено в {out_path} -- проверьте IPRINT (нужно >=2) и что "
              f"расчёт вообще завершился (см. конец файла на предмет "
              f"'*** ERROR'/STOP).", file=sys.stderr)
    return df


# =========================================================================
# 3c. ПАРСИНГ РЕЗУЛЬТАТА ПОДСЧЁТА УЗЛОВ (count-resonances)
# =========================================================================
#
# Формат -- из ИСХОДНИКА fld.driver.f (FORMAT 2410/2420), см. комментарий
# перед NODE_COUNT_TEMPLATE:
#   "  NO NODES BETWEEN FLDMAX AND FLDMIN." -- 0 резонансов в окне.
#   "  NODE COUNT INCREASES BETWEEN FLDMIN AND FLDMAX\n"
#   "  PROGRAM WILL ASSUME MONOTONIC BEHAVIOUR AND SEEK   5 STATES IN THIS
#    INTERVAL" -- 5 резонансов ("INCREASES"/"DECREASES" -- направление
#   счёта узлов с ростом поля, само число резонансов от направления не
#   зависит, это |NODHI-NODLO|).
_RE_NODE_COUNT_CHANGE = re.compile(
    r"NODE COUNT\s+(IN|DE)CREASES BETWEEN FLDMIN AND FLDMAX\s+"
    r"PROGRAM WILL ASSUME MONOTONIC\s+BEHAVIOUR AND SEEK\s+(\d+)\s*STATES?",
    re.IGNORECASE,
)


def parse_node_count(out_path: Path) -> dict:
    """Возвращает {status:'ok', n_resonances: int, direction: 'increases'|
    'decreases'|None} по результату count-resonances (field.exe с MXCALC=3).
    status='no_nodes' (n_resonances=0) для "NO NODES...", status='unknown',
    если ни одно из двух ожидаемых сообщений не нашлось (см. IPRINT/MXCALC
    в NODE_COUNT_TEMPLATE -- возможно, MXCALC оборвал прогон РАНЬШЕ, чем
    успело напечататься сообщение, см. комментарий про порядок операций в
    fld.driver.f)."""
    text = out_path.read_text(errors="replace")
    if _RE_NO_NODES.search(text):
        return {"status": "ok", "n_resonances": 0, "direction": None}
    m = _RE_NODE_COUNT_CHANGE.search(text)
    if m:
        direction = "increases" if m.group(1).upper() == "IN" else "decreases"
        return {"status": "ok", "n_resonances": int(m.group(2)),
                "direction": direction}
    return {"status": "unknown", "n_resonances": None, "direction": None}


def count_resonances_via_field(cfg: Config, fldmin: float, fldmax: float,
                                 monqn: str, label: str) -> dict:
    """Один прогон field.exe (MXCALC=3, см. NODE_COUNT_TEMPLATE) -- БЫСТРЫЙ
    ОТНОСИТЕЛЬНО find_resonances_via_field (2 распространения + максимум 1
    частичное вместо полного перебора N состояний), но НЕ мгновенный:
    стоимость каждого распространения растёт с LMAX (больше каналов в
    сцепленных уравнениях) -- см. переписку: lmax=6 (0.1-20 Гс) уложился в
    12с, lmax=8 -- уже в 76 CPU-с. run_molscat() БЕЗ таймаута -- ждёт
    сколько потребуется, независимо от LMAX/RMAX.
    Возвращает {status, n_resonances, direction, elapsed_sec, out_file}
    БЕЗ поиска точных позиций резонансов."""
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    count_dir = cfg.output_dir / f"{label}_nodecount"
    count_dir.mkdir(parents=True, exist_ok=True)

    print(f"[1/2] Генерация input (MXCALC=3) для field.exe, окно "
          f"{fldmin}-{fldmax} Гс, MONQN={monqn}")
    inp = count_dir / "nodecount.input"
    generate_node_count_input(cfg, fldmin=fldmin, fldmax=fldmax, monqn=monqn,
                               label=f"{label} node count", out_path=inp)

    print("[2/2] Запуск field.exe (только подсчёт узлов на FLDMIN/FLDMAX)...")
    out, elapsed = run_molscat(cfg, inp)
    print(f"       -> заняло {_fmt_duration(elapsed)}")

    result = parse_node_count(out)
    result["elapsed_sec"] = round(elapsed, 1)
    result["out_file"] = str(out)

    if result["status"] == "ok":
        if result["n_resonances"] == 0:
            print(f"       -> 0 резонансов в окне {fldmin}-{fldmax} Гс "
                  f"(NO NODES BETWEEN FLDMAX AND FLDMIN)")
        else:
            print(f"       -> {result['n_resonances']} резонанс(ов) в окне "
                  f"{fldmin}-{fldmax} Гс (число узлов "
                  f"{result['direction']} с ростом поля)")
    else:
        print(f"       -> НЕ УДАЛОСЬ определить число резонансов -- "
              f"ни 'NO NODES...', ни 'NODE COUNT...SEEK N STATES' не "
              f"найдены в {out}. Проверьте .out глазами (возможно, "
              f"MXCALC=3 всё равно оборвал прогон раньше печати -- см. "
              f"комментарий перед NODE_COUNT_TEMPLATE).", file=sys.stderr)
    return result


# =========================================================================
# 4. АВТОДЕТЕКЦИЯ КАНДИДАТОВ
# =========================================================================

def _detect_candidates_in_series(fields: np.ndarray, re_a: np.ndarray, im_a: np.ndarray,
                                   im_a_threshold: float, re_a_jump_sigma: float,
                                   min_separation_G: float) -> list[dict]:
    """Детекция кандидатов внутри ОДНОГО канала (один и тот же физический
    смысл величины re_a/im_a на всём протяжении -- см. detect_candidates)."""
    candidates = []

    abs_im = np.abs(im_a)
    n = len(abs_im)
    # ВАЖНО (найдено на реальных данных, резонанс ~180 Гс пропущен в первом
    # прогоне): канал появляется в таблице "K-DEPENDENT..." ТОЛЬКО когда он
    # достаточно близок к своему порогу (низкая кинетическая энергия) --
    # значит его временной ряд может НАЧИНАТЬСЯ/ЗАКАНЧИВАТЬСЯ ровно на пике
    # резонанса, без соседа с одной из сторон. Раньше цикл шёл по
    # range(1, n-1) и такие краевые точки вообще не проверялись. Теперь
    # проверяем все точки, а отсутствующего соседа считаем "не мешающим"
    # (граничное условие всегда выполнено).
    for i in range(n):
        left_ok = (i == 0) or (abs_im[i] >= abs_im[i - 1])
        right_ok = (i == n - 1) or (abs_im[i] >= abs_im[i + 1])
        if abs_im[i] > im_a_threshold and left_ok and right_ok:
            candidates.append({
                "field_center": fields[i],
                "reason": f"|Im(a)| peak = {abs_im[i]:.2f}"
                          + (" (край ряда канала)" if i == 0 or i == n - 1 else ""),
            })

    if len(re_a) >= 7:
        med = pd.Series(re_a).rolling(7, center=True, min_periods=3).median()
        resid = np.abs(re_a - med.values)
        thr = re_a_jump_sigma * np.nanmedian(resid[resid > 0]) if np.any(resid > 0) else np.inf
        for i in range(len(re_a)):
            if np.isfinite(resid[i]) and resid[i] > thr and thr > 0:
                candidates.append({
                    "field_center": fields[i],
                    "reason": f"Re(a) jump, resid={resid[i]:.2f}",
                })

    candidates.sort(key=lambda c: c["field_center"])
    merged = []
    for c in candidates:
        if merged and abs(c["field_center"] - merged[-1]["field_center"]) < min_separation_G:
            merged[-1]["reason"] += " + " + c["reason"]
            continue
        merged.append(c)
    return merged


def detect_candidates(df: pd.DataFrame, im_a_threshold: float = 0.5,
                       re_a_jump_sigma: float = 5.0,
                       min_separation_G: float = 15.0,
                       l_filter: Optional[int] = 0,
                       window_mult: float = 0.25) -> list[dict]:
    """Резонанс-кандидат = локальный пик |Im(a)| ИЛИ резкая смена знака/скачок
    Re(a) относительно локального фона. Возвращает список окон-кандидатов
    {field_center, field_lo, field_hi, reason, channel, L}.

    ВАЖНО (найдено на реальных данных, см. переписку про 495 Гс/1757 Гс):
    таблица "K-DEPENDENT SCATTERING LENGTHS" печатает РАЗНЫЕ физические
    величины для разных L -- для L=0 это настоящая длина рассеяния (Å,
    POW=1 в calca.f), а для L=2,4 это обобщённые "гиперобъёмы" другой
    размерности (POW=4 и выше) со значениями на порядки больше -- их
    величина НЕ сравнима напрямую с L=0 и почти всегда выглядит как
    "гигантский выброс", даже когда никакого резонанса Фешбаха там нет.
    Поэтому по умолчанию (l_filter=0) детекция кандидатов идёт ТОЛЬКО по
    L=0 каналам (стандартная практика для s-волновых ультрахолодных
    столкновений). l_filter=None включает все L (тогда сравнение амплитуд
    между разными L не имеет физического смысла и легко даёт ложные
    срабатывания -- используй с осторожностью, скорее для отладки).

    Работает на МУЛЬТИКАНАЛЬНОМ df (несколько строк на одно поле, по
    одной на каждый CHAN/L) -- ищет кандидатов НЕЗАВИСИМО в каждом канале
    (одна физическая величина = один временной ряд по полю), затем
    объединяет результаты. Каждый возвращённый кандидат несёт свой
    правильный `channel` (= ICHAN для последующей характеризации) --
    ЭТО ВАЖНО: разные резонансы в разных полевых диапазонах могут сидеть
    в разных фиксированных номерах канала (см. диалог: ICHAN=7 для
    495 Гс, ICHAN=5 для 1757 Гс -- единого on ICHAN на весь диапазон
    0-2000 Гс не существует).

    min_separation_G и window_mult теперь ОБА настраиваемые (см. переписку
    про пилотный скан 0-20 Гс с шагом 0.2 Гс): дефолты 15.0 Гс/0.25
    откалиброваны под грубую сетку 20 Гс и НЕ подходят при мелком шаге --
    min_separation_G=15 на диапазоне 0-20 Гс склеит почти весь скан в
    одного кандидата, а окно field_lo/hi = center +/- dfield/4 при
    dfield=0.2 станет ±0.05 Гс -- слишком узким для надёжной сходимости
    LOCPOL. Явно задавайте оба параметра при отклонении от стандартного
    шага 20 Гс (CLI: --min-separation-g, --window-mult)."""
    if df.empty:
        return []

    if l_filter is not None:
        df = df[df.L == l_filter]
    if df.empty:
        print(f"[WARN] После фильтра L={l_filter} не осталось строк -- "
              f"либо в этом диапазоне поля нет открытых L={l_filter} каналов "
              f"с низкой кинетической энергией, либо l_filter выставлен "
              f"неверно.", file=sys.stderr)
        return []

    all_candidates = []
    for channel, sub in df.groupby("channel"):
        sub = sub.sort_values("field_G").reset_index(drop=True)
        if len(sub) < 3:
            continue  # мало точек для скользящего окна/детекции пика
        found = _detect_candidates_in_series(
            sub.field_G.values, sub.re_a.values, sub.im_a.values,
            im_a_threshold, re_a_jump_sigma, min_separation_G,
        )
        for c in found:
            c["channel"] = int(channel)
            c["L"] = int(sub.L.iloc[0])
            all_candidates.append(c)

    # merge кандидатов из РАЗНЫХ каналов, если поля совпадают в пределах
    # min_separation_G (скорее всего один и тот же физический резонанс,
    # видимый одновременно в нескольких почти вырожденных L=0 каналах) --
    # оставляем ОДНУ запись (первую), не отбрасывая, если поля НЕ совпадают
    # (это могут быть два разных, перекрывающихся резонанса -- см. грабли
    # v3, §4, про 518-536 Гс и 1460-1520 Гс, там как раз нужно НЕ схлопывать).
    all_candidates.sort(key=lambda c: c["field_center"])
    merged = []
    for c in all_candidates:
        if merged and abs(c["field_center"] - merged[-1]["field_center"]) < min_separation_G:
            merged[-1]["reason"] += f" (+ канал {c['channel']}: {c['reason']})"
            continue
        merged.append(c)

    dfield = np.median(np.diff(df.field_G.unique())) if df.field_G.nunique() > 1 else 20.0
    # Узкое окно (см. грабли v3, §3: широкое окно рискует увести экстраполяцию
    # IFCONV в область с другим NOPEN). Пользовательский рабочий пример для
    # кандидата ~495 Гс при шаге сетки 20 Гс использовал окно 490-500
    # (=центр +/- dfield/4, т.е. window_mult=0.25) -- это дефолт, но при
    # мелком шаге сетки (см. пилотный скан 0-20/0.2) окно ±dfield/4
    # становится нереалистично узким -- увеличивайте window_mult (CLI:
    # --window-mult), ориентир -- окно должно быть заметно шире типичной
    # ожидаемой ширины Delta, а не долей шага сетки.
    for c in merged:
        c["field_lo"] = c["field_center"] - dfield * window_mult
        c["field_hi"] = c["field_center"] + dfield * window_mult

    return merged


# =========================================================================
# 5. ХАРАКТЕРИЗАЦИЯ КАНДИДАТА (IFCONV, автоэскалация)
# =========================================================================
#
# Формат подтверждён по source_code/locpol.f (MOLSCAT 2025.0):
#
# ВАЖНО: mol.driver.f вызывает `CALL LOCPOL(...,IFCONV-1,...)`, т.е.
#   IDECAY = IFCONV - 1:
#     IFCONV=1 -> IDECAY=0 (чисто упругий)
#     IFCONV=2 -> IDECAY=1 (слабая фоновая неупругость)
#     IFCONV=3 -> IDECAY=2 (полный неупругий, через DECFIT)
#
# При IPRINT=6 (наш стандартный режим):
#   - IDECAY=0 и IDECAY=1: WRITE(6,130) на КАЖДОЙ итерации даёт B0/DELTA/A_BG
#     (без A_RES/GAMMA_INEL). Полная разбивка с A_RES/GAMMA_INEL для IDECAY=1
#     (FORMAT 170) печатается ТОЛЬКО если IPRINT в диапазоне [3,5] -- при
#     IPRINT=6 её не будет!
#   - IDECAY=2: WRITE(6,106) на КАЖДОЙ итерации даёт ПОЛНЫЙ набор
#     (B0, A_BG re/im, A_RES re/im, GAMMA_INEL, DELTA) БЕЗ ограничения по
#     диапазону IPRINT -- всегда доступно при IPRINT>=6.
#
# ПОЭТОМУ конвейер использует стратегию:
#   IFCONV=1 -> если сходится чисто (не decaying) -- готово (a_res/gamma=None).
#   Если LOCPOL сам сообщает "TERMINATING CONVERGENCE... IFCONV=2 OR 3" ->
#   пропускаем IFCONV=2 и сразу пробуем IFCONV=3 (чтобы получить ПОЛНЫЙ набор
#   параметров одним прогоном, а не рисковать неполными данными на IFCONV=2
#   при IPRINT=6). Если нужно принципиально попробовать IFCONV=2 (например,
#   резонанс на границе применимости DECFIT) -- задай cfg.try_ifconv2=True,
#   тогда конвейер попробует 2 первым, но БЕЗ a_res/gamma (см. выше), и в
#   любом случае эскалирует на 3, если 2 не дал финальной сходимости.

RE_ESCALATE = re.compile(
    r"SHOULD BE CHARACTERISED USING\s+IFCONV\s*=\s*2\s*OR\s*3",
    re.IGNORECASE,
)
RE_NOPEN_CHANGED = re.compile(r"NOPEN CHANGED BETWEEN CALCULATIONS", re.IGNORECASE)
RE_OSCILLATING = re.compile(r"PROBABLY OSCILLATING SO ABANDON CONVERGENCE", re.IGNORECASE)
# "RESONANCE CHARACTERISATION NOT ACHIEVED IN  20 STEPS" -- алгоритм честно
# отработал лимит итераций (по умолчанию 20), но оценка B_RES не
# стабилизировалась (см. переписку про 260/440 Гс: B_RES монотонно уезжал
# на протяжении всех итераций). ЭТО НЕ NOPEN CHANGED и не "oscillating" в
# смысле locpol.f -- отдельная, более тревожная причина: скорее всего это
# НЕ изолированный простой резонанс (наложение нескольких близких резонансов,
# слишком широкий/пологий фон, или иная патология), и трёхточечная формула
# полюса физически не может для него сойтись -- сужение окна тут вряд ли
# поможет, нужен ручной разбор (мелкий локальный скан + визуальный осмотр).
RE_NOT_ACHIEVED = re.compile(
    r"RESONANCE CHARACTERISATION NOT ACHIEVED IN\s*(\d+)\s*STEPS", re.IGNORECASE
)

# "CONVERGED ON RESONANCE AT MAGNETIC Z FIELD_RES = 495.388770 GAUSS, WITH
#  PREDICTED STEP = 5.56622E-04 GAUSS" -- ИМЯ ПЕРЕМЕННОЙ ПЕРЕД "_RES" МОЖЕТ
# СОСТОЯТЬ ИЗ НЕСКОЛЬКИХ СЛОВ ("MAGNETIC Z FIELD"), поэтому используем
# нежадный ".*?" вместо "\S+" (подтверждено реальным выводом пользователя).
RE_CONVERGED = re.compile(
    rf"CONVERGED ON RESONANCE AT\s+.*?_RES\s*=\s*({_NUM})\s*(\S+)"
    rf",\s*WITH PREDICTED STEP\s*=\s*({_NUM})\s*(\S+)",
    re.IGNORECASE,
)

# FORMAT 130 (IDECAY<2): три строки, например
#   "  3-POINT POLE FORMULA ESTIMATES              B_RES = 495.389019328 GAUSS"
#   "                                               DELTA = 7.51000      GAUSS"
#   "                                                A_BG = 47.150       ANGSTROM"
RE_POLE_LINE1 = re.compile(rf"3-POINT POLE FORMULA ESTIMATES.*?=\s*({_NUM})\s*(\S+)", re.IGNORECASE)
RE_POLE_DELTA = re.compile(rf"\bDELTA\s*=\s*({_NUM})\s*(\S+)", re.IGNORECASE)
RE_POLE_ABG = re.compile(rf"\bA_BG\s*=\s*({_NUM})\s*(\S+)", re.IGNORECASE)

# FORMAT 130 (IDECAY=3, т.е. IFCONV=4, Брейт-Вигнер по сумме фаз): три строки,
# ПОДТВЕРЖДЕНО на реальном выводе (440 Гс, 15.07.2026):
#   "  3-POINT POLE FORMULA ESTIMATES MAGNETIC Z FIELD_RES =   444.470966     GAUSS"
#   "                                                GAMMA =  -1.1109         GAUSS"
#   "                                        EPSUM_BG / PI = -0.38329    "
# (третья строка БЕЗ единиц измерения -- EPSUM_BG безразмерна, в долях pi)
RE_POLE_GAMMA = re.compile(rf"\bGAMMA\s*=\s*({_NUM})\s*(\S+)", re.IGNORECASE)
RE_POLE_EPSUM_BG = re.compile(rf"EPSUM_BG\s*/\s*PI\s*=\s*({_NUM})", re.IGNORECASE)
#   "  PARAMETERS OBTAINED ARE:"
#   "  B_RES = 495.389019328 GAUSS"
#   "                             A_BG  (IN ANGSTROM) = 47.150      -1.7500      i "
#   "                             A_RES (IN ANGSTROM) = 214.56       5.8800      i "
#   ""
#   "                                  GAMMA_INEL = -3.3000      GAUSS"
#   "  DELTA (= -ALPHA_RES*GAMMA_INEL/2*ALPHA_BG) = 7.5100      GAUSS"
RE_PARAMS_HEADER = re.compile(r"PARAMETERS OBTAINED ARE:", re.IGNORECASE)
RE_FULL_BRES = re.compile(rf".*?_RES\s*=\s*({_NUM})\s*(\S+)", re.IGNORECASE)
RE_FULL_ABG = re.compile(
    rf"A_BG\s*\(IN\s*(\S+)\)\s*=\s*({_NUM})\s+([+\-]{_NUM})\s*i", re.IGNORECASE
)
RE_FULL_ARES = re.compile(
    rf"A_RES\s*\(IN\s*\S+\)\s*=\s*({_NUM})\s+([+\-]{_NUM})\s*i", re.IGNORECASE
)
RE_FULL_GAMMA = re.compile(rf"GAMMA_INEL\s*=\s*({_NUM})\s*(\S+)", re.IGNORECASE)
RE_FULL_DELTA = re.compile(rf"DELTA\s*\([^)]*\)\s*=\s*({_NUM})\s*(\S+)", re.IGNORECASE)


def _last_match(pattern: re.Pattern, text: str):
    matches = list(pattern.finditer(text))
    return matches[-1] if matches else None


def parse_ifconv_block(out_path: Path, ifconv: int) -> dict:
    """Парсит .out одной характеризации. Возвращает:
    {status: 'ok'|'escalate'|'nopen_changed'|'oscillating'|'unknown',
     B0_G, Delta_G, a_bg_re, a_bg_im, a_res_re, a_res_im, Gamma_inel_G
     (последние 4 присутствуют только для ifconv==3, см. комментарий выше)}"""
    text = out_path.read_text(errors="replace")

    if RE_NOPEN_CHANGED.search(text):
        return {"status": "nopen_changed"}
    if RE_OSCILLATING.search(text):
        return {"status": "oscillating"}
    if RE_NOT_ACHIEVED.search(text):
        result = {"status": "not_converged"}
        # НЕ сошлось строго по DTOL, но последняя оценка часто уже полезна
        # (см. переписку про 440 Гс IFCONV=4: B_RES колебался в узкой полосе
        # 441-444 Гс, несмотря на формальную "неудачу") -- сохраняем её как
        # tentative, явно помечая, что это не финальный результат.
        m_pole = _last_match(RE_POLE_LINE1, text)
        if m_pole:
            result["tentative_B0_G"] = float(m_pole.group(1))
            tail = text[m_pole.end():m_pole.end() + 300]
            if ifconv == 4:
                m_gamma = RE_POLE_GAMMA.search(tail)
                if m_gamma:
                    result["tentative_Gamma_BW_G"] = float(m_gamma.group(1))
            else:
                m_delta = RE_POLE_DELTA.search(tail)
                if m_delta:
                    result["tentative_Delta_G"] = float(m_delta.group(1))
        return result
    if ifconv == 1 and RE_ESCALATE.search(text):
        return {"status": "escalate"}

    m_conv = RE_CONVERGED.search(text)
    if not m_conv:
        return {"status": "unknown"}

    result = {"status": "ok", "B0_G": float(m_conv.group(1)),
              "predicted_step": float(m_conv.group(3))}

    if ifconv == 3:
        # ищем ПОСЛЕДНИЙ блок "PARAMETERS OBTAINED ARE:" (финальная итерация)
        m_hdr = _last_match(RE_PARAMS_HEADER, text)
        if m_hdr:
            tail = text[m_hdr.end():m_hdr.end() + 1000]
            m_bres = RE_FULL_BRES.search(tail)
            m_abg = RE_FULL_ABG.search(tail)
            m_ares = RE_FULL_ARES.search(tail)
            m_gamma = RE_FULL_GAMMA.search(tail)
            m_delta = RE_FULL_DELTA.search(tail)
            if m_bres:
                result["B0_G"] = float(m_bres.group(1))
            if m_abg:
                result["a_bg_re"] = float(m_abg.group(2))
                result["a_bg_im"] = float(m_abg.group(3))
            if m_ares:
                result["a_res_re"] = float(m_ares.group(1))
                result["a_res_im"] = float(m_ares.group(2))
            if m_gamma:
                result["Gamma_inel_G"] = float(m_gamma.group(1))
            if m_delta:
                result["Delta_G"] = float(m_delta.group(1))
    elif ifconv == 4:
        # IDECAY=3, Брейт-Вигнер по сумме собственных фаз: та же строка-якорь
        # "3-POINT POLE FORMULA ESTIMATES ..._RES = ...", но следом идут
        # GAMMA (ширина Брейта-Вигнера по ПОЛНОЙ сумме собственных фаз --
        # это ровно Gamma_B из Hutson, NJP 9, 152 (2007), Eq. (6): для
        # ОДНОГО открытого канала сумма фаз -- это просто сам фазовый сдвиг
        # delta_0, так что это тот же формализм, что и Eq. (5)-(9) статьи,
        # НЕ равна Delta_G напрямую без пересчёта) и безразмерная EPSUM_BG
        # (фоновая сумма фаз / pi, т.е. delta_bg/pi).
        m_pole = _last_match(RE_POLE_LINE1, text)
        if m_pole:
            tail = text[m_pole.end():m_pole.end() + 300]
            m_gamma = RE_POLE_GAMMA.search(tail)
            m_epsbg = RE_POLE_EPSUM_BG.search(tail)
            if m_gamma:
                result["Gamma_BW_G"] = float(m_gamma.group(1))
            if m_epsbg:
                result["EPSUM_bg_over_pi"] = float(m_epsbg.group(1))
        result["a_res_re"] = None
        result["Gamma_inel_G"] = None

        # ПЕРЕСЧЁТ Delta_B из Gamma_B: Hutson NJP 9, 152 (2007), Eq. (2) и
        # (9) -- a_bg = -tan(delta_bg)/k, Gamma_B = -2*a_bg*k*Delta_B.
        # Подставляя delta_bg = pi*EPSUM_bg_over_pi (одноканальный случай),
        # a_bg*k = -tan(pi*EPSUM_bg_over_pi), и k СОКРАЩАЕТСЯ:
        #   Delta_B = Gamma_B / (2*tan(pi*EPSUM_bg_over_pi))
        # Проверено эмпирически в toy_resonance_pipeline.py (state 8.93 Гс,
        # lmax=10): форсированный IFCONV=4 на резонансе с уже известным
        # Delta_G=0.0112959 (из IFCONV=1) дал Delta_B(пересчёт)=0.011526 --
        # согласие ~2%. ВАЖНО: формула строго верна только для ОДНОГО
        # открытого канала (иначе EPSUM_bg_over_pi -- сумма ФАЗ нескольких
        # каналов, а не delta_bg одного канала).
        if result.get("Gamma_BW_G") is not None and result.get("EPSUM_bg_over_pi") is not None:
            tan_bg = math.tan(math.pi * result["EPSUM_bg_over_pi"])
            if tan_bg != 0:
                result["Delta_G_from_gamma_bw"] = result["Gamma_BW_G"] / (2 * tan_bg)
    else:
        # IFCONV 1 или 2: FORMAT 130 -- B0/DELTA/A_BG, без a_res/gamma
        m_pole = _last_match(RE_POLE_LINE1, text)
        if m_pole:
            tail = text[m_pole.end():m_pole.end() + 300]
            m_delta = RE_POLE_DELTA.search(tail)
            m_abg = RE_POLE_ABG.search(tail)
            if m_delta:
                result["Delta_G"] = float(m_delta.group(1))
            if m_abg:
                result["a_bg_re"] = float(m_abg.group(1))
        result["a_res_re"] = None
        result["Gamma_inel_G"] = None

    # ГАРАНТИЯ КОНТРАКТА (см. докстринг выше): вызывающий код
    # (characterize_candidate / build_summary_table) ожидает, что при
    # status=="ok" присутствуют ВСЕ перечисленные ключи, хотя реально их
    # печатает LOCPOL только при определённых условиях -- например, для
    # ЧИСТО УПРУГОГО резонанса, сошедшегося уже на IFCONV=1 (случай "если
    # сходится чисто -- готово" из докстринга characterize_candidate),
    # блок "3-POINT POLE FORMULA ESTIMATES ... DELTA = ..." в .out вообще
    # не печатается (неупругих величин там просто нет), так что m_pole/
    # m_delta выше не находятся и "Delta_G" не проставляется. Раньше это
    # приводило к KeyError в build_summary_table (r["Delta_G"]) на любом
    # чисто упругом резонансе. Проставляем None всем ожидаемым, но ещё не
    # заданным полям, чтобы результат всегда был однородным по набору
    # ключей независимо от того, какие строки реально нашлись в .out.
    for _key in ("Delta_G", "a_bg_re", "a_bg_im", "a_res_re", "a_res_im",
                 "Gamma_inel_G", "Gamma_BW_G", "EPSUM_bg_over_pi",
                 "Delta_G_from_gamma_bw"):
        result.setdefault(_key, None)

    return result


def characterize_candidate(cfg: Config, candidate: dict, run_id: str) -> dict:
    """Полный цикл характеризации кандидата.

    Стратегия (см. комментарий над регэкспами выше про асимметрию IPRINT
    между IFCONV=2 и IFCONV=3 в реальном коде locpol.f):
      1. IFCONV=1. Если сходится чисто упруго -- готово.
      2. Если LOCPOL сам просит эскалацию явным сообщением -- пробуем СРАЗУ
         IFCONV=3 (даёт полный набор параметров одним прогоном при IPRINT=6).
      3. NOPEN CHANGED / OSCILLATING / NOT_CONVERGED -- сужаем бракетирующее
         окно вдвое (до 3 раз) и повторяем ТОТ ЖЕ уровень IFCONV.
      4. Если сужения на текущем уровне исчерпаны, А ОН ТАК И НЕ СОШЁЛСЯ --
         эскалируем на следующий уровень (1->3->4) на ИСХОДНОМ окне
         кандидата, даже если LOCPOL не печатал явную просьбу об эскалации
         (см. ниже -- это ИСПРАВЛЕНИЕ реального бага, портированное из
         toy_resonance_pipeline.py: раньше эскалация 1->3 срабатывала
         ТОЛЬКО по явному сообщению LOCPOL, и если на IFCONV=1 после 3
         сужений статус был not_converged/nopen_changed/oscillating, а не
         'escalate', конвейер сдавался после ровно 4 вызовов molscat, ни
         разу не попробовав IFCONV=3/4). IFCONV=4 (Брейт-Вигнер по сумме
         собственных фаз) -- принципиально другая, более устойчивая
         величина (см. переписку про 440 Гс), может сойтись там, где длина
         рассеяния не смогла. ВАЖНО: лимит MXLOC=20 итераций общий для всех
         IFCONV (жёстко зашит в sizes_module.f, не настраивается через
         &INPUT) -- смена уровня не даёт "больше попыток" САМОМУ LOCPOL,
         только другую точку входа/величину для той же 3-точечной
         экстраполяции.
      5. Если и IFCONV=4 не сошёлся после всех сужений -- возвращаем
         'failed', но сохраняем ЛУЧШУЮ предварительную (не до конца
         сошедшуюся) оценку B0 из последней итерации любого из уровней --
         это по-прежнему полезная информация (см. 440 Гс: B_RES не сошёлся
         строго, но стабильно лежал в полосе 441-444 Гс). Максимум попыток
         теперь до 12 (4 на каждый из 3 уровней) вместо прежних де-факто 4
         в худшем случае.

    Все .input/.out этого кандидата складываются в ОТДЕЛЬНУЮ подпапку
    cfg.output_dir/<run_id>/ (а не вперемешку с другими кандидатами) --
    у одного резонанса может быть до 12 файлов (эскалации + сужения окна),
    подпапка на кандидата держит это читаемым."""
    cand_dir = cfg.output_dir / run_id
    cand_dir.mkdir(parents=True, exist_ok=True)

    fmin, fmax = candidate["field_lo"], candidate["field_hi"]
    ichan = candidate["channel"]  # каждый кандидат несёт свой правильный ICHAN
                                   # (см. detect_candidates -- разные резонансы
                                   # могут сидеть в разных фиксированных каналах)
    ifconv = 1
    narrow_attempts = 0
    tried_3 = False
    tried_4 = False
    best_tentative = None
    total_elapsed = 0.0
    n_calls = 0

    while True:
        label = f"ifconv{ifconv}_attempt{n_calls+1}"
        inp = cand_dir / f"{label}.input"
        generate_ifconv_input(
            cfg, fmin=fmin, fmax=fmax,
            ichan=ichan, ifconv=ifconv, label=label, out_path=inp,
        )
        out, elapsed = run_molscat(cfg, inp)
        total_elapsed += elapsed
        n_calls += 1
        print(f"          [{_fmt_duration(elapsed)}] IFCONV={ifconv}, "
              f"окно {fmin:.2f}-{fmax:.2f} Гс -> {out.name}")
        result = parse_ifconv_block(out, ifconv)

        if result["status"] == "ok":
            result["candidate"] = candidate
            result["ifconv_used"] = ifconv
            result["out_file"] = str(out)
            result["elapsed_sec"] = round(total_elapsed, 1)
            result["n_molscat_calls"] = n_calls
            return result

        # Явная просьба LOCPOL (только на IFCONV=1) -- эскалируем СРАЗУ,
        # без исчерпания сужений на IFCONV=1.
        if result["status"] == "escalate" and not tried_3:
            tried_3 = True
            ifconv = 3
            narrow_attempts = 0
            continue

        if result.get("tentative_B0_G") is not None:
            best_tentative = {**result, "ifconv_used": ifconv}

        if result["status"] in ("nopen_changed", "oscillating", "not_converged") and narrow_attempts < 3:
            narrow_attempts += 1
            span = (fmax - fmin) / 2
            fmin = candidate["field_center"] - span / 2
            fmax = candidate["field_center"] + span / 2
            continue

        # Сужения на ТЕКУЩЕМ уровне исчерпаны (или статус вообще не входит
        # в список выше, напр. 'unknown') -- эскалируем на следующий
        # уровень 1->3->4 на исходном окне кандидата, ВМЕСТО немедленной
        # сдачи (см. докстринг выше -- это исправление реального бага).
        if not tried_3:
            tried_3 = True
            ifconv = 3
            fmin, fmax = candidate["field_lo"], candidate["field_hi"]
            narrow_attempts = 0
            continue
        if not tried_4:
            tried_4 = True
            ifconv = 4
            fmin, fmax = candidate["field_lo"], candidate["field_hi"]
            narrow_attempts = 0
            continue

        failed = {"status": "failed", "last_status": result["status"],
                  "candidate": candidate, "out_file": str(out),
                  "elapsed_sec": round(total_elapsed, 1), "n_molscat_calls": n_calls}
        if best_tentative is not None:
            failed["tentative"] = best_tentative
        return failed


# =========================================================================
# 5b. ПОИСК ICHAN ПО (F,mF): THRSH9 сам НЕ нумерует каналы -- он только
# считает число EREF (см⁻¹) по Брейту-Раби для заданных MONQN. Номер
# канала (ICHAN) -- это позиция строки в таблице "THRESHOLDS CALCULATED
# FROM ASYMPTOTIC HAMILTONIAN" конкретного прогона (JTOT/IBFIX/поле),
# которую нужно найти по СОВПАДЕНИЮ ЧИСЕЛ с напечатанным EREF. Раньше
# это делалось вручную (см. переписку про ICHAN=7 vs 5 для разных
# резонансов); ниже -- автоматизация того же самого сопоставления.
# =========================================================================

# "  REFERENCE ENERGY IS                          -1.5238287556E-02 CM-1    = ..."
RE_REFERENCE_ENERGY = re.compile(
    rf"REFERENCE ENERGY IS\s+({_NUM})\s*CM-1", re.IGNORECASE
)

# "  THRESHOLD       L           ENERGY/CM-1                ENERGY/K"
RE_THRESHOLD_HEADER = re.compile(r"THRESHOLD\s+L\s+ENERGY/CM-1", re.IGNORECASE)

# "          1       0        -1.523828755595E-02        -2.192449579097E-02"
RE_THRESHOLD_ROW = re.compile(
    rf"^\s*(\d+)\s+(\d+)\s+({_NUM})\s+({_NUM})\s*$", re.MULTILINE
)


def parse_threshold_table(out_path: Path) -> dict:
    """Разбирает ПЕРВЫЙ встреченный блок 'REFERENCE ENERGY IS ...' +
    'THRESHOLDS CALCULATED FROM ASYMPTOTIC HAMILTONIAN' в .out (для
    диагностического прогона find_channel_for_monqn их всегда ровно
    один, т.к. FLDMIN=FLDMAX -- одна точка по полю).
    Возвращает {reference_energy_cm1, thresholds: [{index,L,energy_cm1,
    energy_K}, ...]} или {} если не нашлось (см. подсказки в
    find_channel_for_monqn про IPRINT<6)."""
    text = out_path.read_text(errors="replace")

    m_ref = RE_REFERENCE_ENERGY.search(text)
    if not m_ref:
        return {}
    reference_energy_cm1 = float(m_ref.group(1).replace("D", "E").replace("d", "e"))

    m_hdr = RE_THRESHOLD_HEADER.search(text, pos=m_ref.end())
    if not m_hdr:
        return {"reference_energy_cm1": reference_energy_cm1, "thresholds": []}

    # таблица заканчивается на первой пустой строке после заголовка
    tail = text[m_hdr.end():]
    end = tail.find("\n\n")
    block = tail if end < 0 else tail[:end]

    thresholds = []
    for row in RE_THRESHOLD_ROW.finditer(block):
        thresholds.append({
            "index": int(row.group(1)),
            "L": int(row.group(2)),
            "energy_cm1": float(row.group(3).replace("D", "E").replace("d", "e")),
            "energy_K": float(row.group(4).replace("D", "E").replace("d", "e")),
        })

    return {"reference_energy_cm1": reference_energy_cm1, "thresholds": thresholds}


def find_channel_for_monqn(cfg: Config, monqn: str, field: float,
                            label: str, tol_cm1: float = 1.0e-8,
                            l_filter: Optional[int] = 0,
                            out_dir: Optional[Path] = None) -> dict:
    """Гоняет ОДНОТОЧЕЧНЫЙ диагностический прогон (FLDMIN=FLDMAX=field) с
    заданным MONQN (через cfg.monqn) и находит номер канала (строку в
    таблице THRESHOLDS), чья энергия совпадает с EREF, посчитанным
    THRSH9 по этим F,mF. ВАЖНО: как обсуждали -- при B=0 разные (F,mF)
    комбинации могут быть точно вырождены (несколько строк с одинаковой
    энергией), поэтому по умолчанию используйте НЕНУЛЕВОЕ поле (см.
    field).

    ВАЖНО #2 (портировано из toy_resonance_pipeline.py -- ambiguous при
    LMAX>=2 на ЛЮБОМ поле, не только вблизи вырождения уровней): энергия
    порога в THRSH9/THRESHOLDS зависит ТОЛЬКО от внутреннего (F,mF)-
    состояния, не от партциальной волны L -- значит для одного и того же
    MONQN таблица THRESHOLDS всегда содержит ПО ОДНОЙ строке на каждый
    открытый L (L=0,2,4,...) с ОДИНАКОВОЙ energy_cm1, даже вдали от любых
    реальных пересечений уровней. l_filter (дефолт 0, s-волна -- см. тот
    же дефолт в detect_candidates/l_filter) разрешает эту неоднозначность:
    если среди совпавших по энергии строк ровно одна имеет L==l_filter,
    она возвращается как единственный результат (без статуса 'ambiguous').
    Если l_filter=None или после фильтра всё ещё >1 строки -- ведёт себя
    как раньше (см. ниже), сохраняя ВСЕ совпавшие строки в 'candidates'
    для ручной диагностики (это уже настоящая неоднозначность, обычно
    из-за близости к точке пересечения уровней).

    out_dir (по умолчанию cfg.output_dir) -- куда писать .input/.out;
    вызывающий код (characterize_field_resonances) передаёт СВОЮ подпапку
    на резонанс, чтобы *_chan.input/.out не мусорили верхний уровень
    output_dir (см. переписку про toy_1_50000G_4/6 -- сотни файлов
    *_chan.* напрямую в molscat_runs).

    Возвращает:
      status='ok': {channel, L, energy_cm1, diff_cm1, ...}
      status='ambiguous': несколько строк совпали в пределах tol_cm1 (после
        применения l_filter, если он задан) -- обычно значит, что field
        слишком близко к точке вырождения (НЕ обычная L-деградация, см.
        выше -- та уже разрешена); candidates перечисляет все совпавшие
        строки.
      status='not_found': ни одна строка не совпала -- см. all_thresholds
        для ручной диагностики (например, если MONQN задан некорректно,
        или IPRINT<6 в шаблоне, или ANSA=0 (INUCA=0) -- любой из этих
        случаев даст пустой/бессмысленный список)."""
    if not monqn:
        raise ValueError("find_channel_for_monqn: monqn не задан")
    cfg_local = dataclasses.replace(cfg, monqn=monqn)
    target_dir = out_dir if out_dir is not None else cfg.output_dir
    target_dir.mkdir(parents=True, exist_ok=True)
    inp = target_dir / f"{label}.input"
    generate_coarse_scan_input(
        cfg_local, fmin=field, fmax=field, dfield=1.0, label=label, out_path=inp,
    )
    out, elapsed = run_molscat(cfg_local, inp)
    parsed = parse_threshold_table(out)

    if not parsed or "reference_energy_cm1" not in parsed:
        return {"status": "not_found", "reason": "REFERENCE ENERGY не найдена в .out "
                "(проверьте, что MONQN синтаксически верен и что THRSH9 "
                "реализован в собранном .exe)", "out_file": str(out)}

    eref = parsed["reference_energy_cm1"]
    thresholds = parsed["thresholds"]
    matches = [t for t in thresholds if abs(t["energy_cm1"] - eref) < tol_cm1]

    if len(matches) > 1 and l_filter is not None:
        l_matches = [t for t in matches if t["L"] == l_filter]
        if len(l_matches) == 1:
            matches = l_matches  # неоднозначность разрешена фильтром по L

    if len(matches) == 1:
        t = matches[0]
        return {"status": "ok", "channel": t["index"], "L": t["L"],
                "energy_cm1": t["energy_cm1"],
                "diff_cm1": t["energy_cm1"] - eref,
                "reference_energy_cm1": eref, "field_G": field,
                "elapsed_sec": round(elapsed, 1), "out_file": str(out)}
    if len(matches) > 1:
        return {"status": "ambiguous", "candidates": matches,
                "reference_energy_cm1": eref, "field_G": field,
                "reason": f"{len(matches)} строк с той же энергией "
                f"(после фильтра L={l_filter}) -- похоже на вырождение "
                "уровней; попробуйте другое поле.", "out_file": str(out)}
    return {"status": "not_found", "all_thresholds": thresholds,
            "reference_energy_cm1": eref, "field_G": field,
            "reason": "ни одна строка THRESHOLDS не совпала с EREF в "
            f"пределах tol_cm1={tol_cm1:.1e} -- см. all_thresholds для "
            "ручного сопоставления (возможно, tol_cm1 слишком мал/велик, "
            "либо MONQN указывает на состояние вне текущего JTOT/IBFIX "
            "сектора).", "out_file": str(out)}


# =========================================================================
# 6. СБОРКА ИТОГОВОЙ ТАБЛИЦЫ
# =========================================================================

def build_summary_table(results: list[dict], meta: RunMetadata) -> pd.DataFrame:
    rows = []
    for r in results:
        if r.get("status") != "ok":
            row = {
                "B0_G": None, "Delta_G": None, "status": r.get("status", "failed"),
                "last_status": r.get("last_status"),
                "field_center_guess": r["candidate"]["field_center"],
                "channel": r["candidate"].get("channel"),
                "L": r["candidate"].get("L"),
                "reason": r["candidate"]["reason"],
                "elapsed_sec": r.get("elapsed_sec"),
                "n_molscat_calls": r.get("n_molscat_calls"),
                "tensor_terms": meta.tensor_terms,
                "hyperfine_zeeman": meta.hyperfine_zeeman,
                "lmax": meta.lmax,
                "dipdip_scale": meta.dipdip_scale,
                "v2_scale": meta.v2_scale,
                "run_label": meta.label,
            }
            tent = r.get("tentative")
            if tent:
                row["tentative_ifconv"] = tent.get("ifconv_used")
                row["tentative_B0_G"] = tent.get("tentative_B0_G")
                row["tentative_Delta_G"] = tent.get("tentative_Delta_G")
                row["tentative_Gamma_BW_G"] = tent.get("tentative_Gamma_BW_G")
            rows.append(row)
            continue
        rows.append({
            "B0_G": r.get("B0_G"), "Delta_G": r.get("Delta_G"),
            "a_bg_re": r.get("a_bg_re"), "a_bg_im": r.get("a_bg_im"),
            "a_res_re": r.get("a_res_re"), "a_res_im": r.get("a_res_im"),
            "Gamma_inel_G": r.get("Gamma_inel_G"),
            "Gamma_BW_G": r.get("Gamma_BW_G"),
            "EPSUM_bg_over_pi": r.get("EPSUM_bg_over_pi"),
            "Delta_G_from_gamma_bw": r.get("Delta_G_from_gamma_bw"),
            "ifconv_used": r.get("ifconv_used"),
            "predicted_step_G": r.get("predicted_step"),
            "channel": r["candidate"].get("channel"),
            "L": r["candidate"].get("L"),
            "field_center_guess": r["candidate"].get("field_center"),
            "elapsed_sec": r.get("elapsed_sec"),
            "n_molscat_calls": r.get("n_molscat_calls"),
            "status": "ok",
            "tensor_terms": meta.tensor_terms,
            "hyperfine_zeeman": meta.hyperfine_zeeman,
            "lmax": meta.lmax,
            "dipdip_scale": meta.dipdip_scale,
            "v2_scale": meta.v2_scale,
            "run_label": meta.label,
        })
    return pd.DataFrame(rows)


def dedupe_by_B0(df: pd.DataFrame, tol_G: float = 2.0) -> pd.DataFrame:
    """Схлопывает дубликаты: РАЗНЫЕ кандидаты грубого скана (обычно с двух
    соседних узлов сетки по разные стороны от резонанса, см. переписку про
    1085.97/1336.38/1757.21 Гс -- LOCPOL стартует с бракетирующего окна
    вокруг узла сетки, но экстраполяция 3-точечной формулы полюса может
    "уехать" на много Гс к настоящему полюсу, если резонанс сидит МЕЖДУ
    двумя соседними точками грубого скана, и оба соседних узла тогда
    независимо сходятся на ОДИН И ТОТ ЖЕ настоящий резонанс). Оставляет
    запись с наименьшим |predicted_step_G| (точнее сошедшуюся), если он
    известен, иначе первую по порядку.

    НЕ путать с реальными близко расположенными, но РАЗНЫМИ резонансами
    (например, 1483.79 и 1515.29 Гс отличаются на ~31 Гс -- намного больше
    tol_G=2 Гс по умолчанию, поэтому останутся отдельными записями)."""
    ok = df[df.status == "ok"].copy()
    other = df[df.status != "ok"]
    if ok.empty:
        return df

    ok = ok.sort_values("B0_G").reset_index(drop=True)
    keep = []
    used = set()
    for i in range(len(ok)):
        if i in used:
            continue
        group = [i]
        for j in range(i + 1, len(ok)):
            if j in used:
                continue
            if abs(ok.loc[j, "B0_G"] - ok.loc[group[-1], "B0_G"]) <= tol_G:
                group.append(j)
            else:
                break
        for idx in group:
            used.add(idx)
        if len(group) == 1:
            keep.append(group[0])
        else:
            # выбираем с наименьшим |predicted_step_G|, если есть данные
            steps = ok.loc[group, "predicted_step_G"]
            if steps.notna().any():
                best = ok.loc[group].loc[steps.abs().idxmin()].name
            else:
                best = group[0]
            keep.append(best)

    return pd.concat([ok.loc[keep], other], ignore_index=True).sort_values(
        "B0_G", na_position="last"
    ).reset_index(drop=True)


# =========================================================================
# 6b. ПОИСК РЕЗОНАНСОВ ЧЕРЕЗ FIELD (альтернатива Этапам 1-5 для
#     одноканальной ситуации -- см. FIELD_TEMPLATE выше)
# =========================================================================

def find_resonances_via_field(cfg: Config, fldmin: float, fldmax: float,
                                monqn: str, label: str,
                                csv_out: Path) -> pd.DataFrame:
    """Один прогон field.exe: ищет ВСЕ поля в [fldmin, fldmax], где связанный
    уровень пересекает порог (ENERGY=0 по MONQN) -- т.е. позиции резонансов
    B0 НАПРЯМУЮ, без грубого скана/detect_candidates/IFCONV. cfg.molscat_exe
    должен указывать на field-Tm2.exe (НЕ molscat-Tm2.exe -- это другая
    программа, см. Makefile: BASE9-TM-TM используется и там, и там, но
    драйверы разные: fld.driver.f против mol.driver.f)."""
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    field_dir = cfg.output_dir / f"{label}_field"
    field_dir.mkdir(parents=True, exist_ok=True)

    print(f"[1/3] Генерация input для field.exe, окно {fldmin}-{fldmax} Гс, "
          f"MONQN={monqn}")
    field_input = field_dir / "field.input"
    generate_field_scan_input(cfg, fldmin=fldmin, fldmax=fldmax, monqn=monqn,
                               label=f"{label} field search", out_path=field_input)

    print("[2/3] Запуск field.exe...")
    field_out, elapsed = run_molscat(cfg, field_input)
    print(f"       -> заняло {_fmt_duration(elapsed)}")

    print("[3/3] Парсинг результатов...")
    df = parse_field_resonances(field_out)
    df.to_csv(csv_out, index=False)

    n_ok = (df.status == "ok").sum() if not df.empty else 0
    n_other = len(df) - n_ok
    if df.empty:
        print(f"       -> резонансов в окне {fldmin}-{fldmax} Гс нет "
              f"(NO NODES BETWEEN FLDMAX AND FLDMIN -- т.е. число связанных "
              f"состояний ниже порога ОДИНАКОВО на обеих границах окна)")
    else:
        print(f"       -> найдено {len(df)} записей ({n_ok} чистых, "
              f"{n_other} с оговоркой -- см. колонку status)")
        for _, r in df.iterrows():
            if r.status == "not_converged":
                print(f"          state {r.state}: НЕ СОШЛОСЬ")
            else:
                print(f"          state {r.state}: B0 = {r.field_G:.6f} "
                      f"{r.unit}  [{r.status}]")
    print(f"       -> записано в {csv_out}")
    return df


# =========================================================================
# 6c. МОСТ FIELD -> LOCPOL/IFCONV (эластичная процедура Frye & Hutson)
# =========================================================================
#
# find_resonances_via_field() даёт точные B0 (позиции резонансов), но САМ
# field.exe не считает ни ширину, ни a_bg -- это отдельная задача (энергия
# связанного состояния относительно порога, не длина рассеяния). Frye & Hutson,
# PRA 96, 042705 (2017), Sec. II ("elastic procedure") описывают ИМЕННО эту
# ситуацию: FIELD даёт стартовое приближение B0 ("we usually use the program
# FIELD..."), а Δ (ширина по полю) и a_bg получаются из 3-точечной
# экстраполяции полюса длины рассеяния, Eq. (2)-(5). Мы это НЕ реализуем
# заново в Python -- она уже есть в locpol.f (см. RE_POLE_LINE1/RE_POLE_DELTA/
# RE_POLE_ABG выше, FORMAT 130 -- буквально B_RES/DELTA/A_BG из статьи) и уже
# вызывается через characterize_candidate() (IFCONV=1 = IDECAY=0 = чисто
# упругий случай статьи; эскалация на IFCONV=3/4 -- их же "regularized"/
# "fully complex" процедуры для случаев с неупругостью, см. комментарии над
# characterize_candidate). Здесь только МОСТ: B0 из field.exe -> ICHAN (через
# find_channel_for_monqn) -> candidate-словарь -> characterize_candidate().
# Портировано из toy_resonance_pipeline.py, вместе со всеми доработками,
# найденными при её обкатке (l_filter в find_channel_for_monqn, подпапка на
# резонанс, эскалация в characterize_candidate, проверка field_match,
# пересчёт Delta_G_from_gamma_bw).

def characterize_field_resonances(cfg_scatter: Config, field_df: pd.DataFrame,
                                    monqn: str, bracket_g: float,
                                    run_prefix: str,
                                    l_filter: Optional[int] = 0) -> list[dict]:
    """Для каждой строки field_df со status=='ok' (надёжный B0 от field.exe):
      1. Определяет ICHAN на этом B0 через find_channel_for_monqn (THRSH9
         сам номер канала не знает, только EREF -- см. секцию 5b).
      2. Собирает candidate-словарь с УЗКИМ бракетирующим окном
         B0 +/- bracket_g (в отличие от кандидатов грубого скана, B0 здесь
         уже точен -- окно нужно только чтобы дать LOCPOL стартовые точки,
         см. Sec. II статьи про δB=0.2 Гс между начальными точками; ЕСЛИ
         окно всё же не подойдёт -- characterize_candidate САМ его сузит,
         см. narrow_attempts в этой функции).
      3. Прогоняет через characterize_candidate() (тот же IFCONV-конвейер
         с эскалацией 1->3->4, что и для кандидатов грубого скана).

    cfg_scatter.molscat_exe ДОЛЖЕН быть molscat-Tm2.exe (характеризация --
    это обычный рассеятельный расчёт с IFCONV, параметром mol.driver.f, а
    не fld.driver.f). Строки со status != 'ok' (not_converged/out_of_range)
    пропускаются с предупреждением -- их B0 недостаточно надёжен, чтобы
    от него стартовать 3-точечную экстраполюцию.

    Возвращает список результатов в формате characterize_candidate() (плюс
    ключи 'state' и 'field_G_field_driver' для прослеживаемости), включая
    отдельный статус 'channel_not_found', если find_channel_for_monqn не
    смог сопоставить канал (само по себе диагностически полезно -- скорее
    всего JTOT/IBFIX сектор не тот, где искал field.exe)."""
    results: list[dict] = []
    n_total = len(field_df)
    for i, row in enumerate(field_df.itertuples(index=False)):
        if row.status != "ok":
            print(f"[WARN] state {row.state}: status={row.status!r} -- "
                  f"пропускаю характеризацию (ненадёжный B0)", file=sys.stderr)
            continue

        b0_guess = row.field_G
        # Общая подпапка на этот резонанс (та же, что чуть ниже использует
        # characterize_candidate) -- сюда же кладём файлы поиска ICHAN,
        # чтобы верхний уровень output_dir не засорялся сотнями *_chan.*.
        run_id = f"{run_prefix}_state{row.state}"
        state_dir = cfg_scatter.output_dir / run_id
        print(f"       -> резонанс {i+1}/{n_total}: state {row.state} @ "
              f"{b0_guess:.6f} Гс -- ищу ICHAN...")
        chan_res = find_channel_for_monqn(cfg_scatter, monqn, b0_guess, "chan",
                                          l_filter=l_filter, out_dir=state_dir)

        if chan_res["status"] != "ok":
            print(f"          ICHAN не определён ({chan_res['status']}: "
                  f"{chan_res.get('reason')}) -- характеризация пропущена",
                  file=sys.stderr)
            results.append({
                "status": "channel_not_found",
                "candidate": {
                    "field_center": b0_guess, "channel": None, "L": None,
                    "reason": f"field.exe state {row.state}: ICHAN lookup "
                              f"failed ({chan_res['status']})",
                },
                "out_file": chan_res.get("out_file"),
                "elapsed_sec": chan_res.get("elapsed_sec"),
                "n_molscat_calls": 1,
                "state": row.state, "field_G_field_driver": b0_guess,
            })
            continue

        candidate = {
            "field_center": b0_guess,
            "field_lo": b0_guess - bracket_g, "field_hi": b0_guess + bracket_g,
            "channel": chan_res["channel"], "L": chan_res["L"],
            "reason": f"field.exe state {row.state}, B0={b0_guess:.6f} Гс",
        }
        print(f"          ICHAN={chan_res['channel']} (L={chan_res['L']}) -- "
              f"характеризация (IFCONV/LOCPOL)...")
        res = characterize_candidate(cfg_scatter, candidate, run_id)
        res["state"] = row.state
        res["field_G_field_driver"] = b0_guess
        if res.get("status") == "ok":
            # ВАЖНО (портировано из toy_resonance_pipeline.py, найдено на
            # реальном прогоне): FLDMIN/FLDMAX в IFCONV -- это ТОЛЬКО
            # стартовые точки 3-точечной экстраполяции (Eq. 2-5 статьи), НЕ
            # жёсткая граница поиска -- если рядом с b0_guess нет НАСТОЯЩЕГО
            # полюса длины рассеяния (напр. этот уровень слабо/не связан с
            # открытым каналом -- field.exe находит ЛЮБОЕ пересечение
            # связанного уровня с порогом, независимо от силы связи, а
            # LOCPOL видит только то, что реально проявляется как
            # особенность a(B)), экстраполяция может "убежать" сколь угодно
            # далеко (видели даже отрицательное поле как промежуточную
            # оценку) и сойтись на СОВСЕМ ДРУГОМ, не связанном с этим state,
            # резонансе. Проверяем это явно, вместо того чтобы молча
            # доверять B0_G/Delta_G/a_bg -- иначе такая ошибочная
            # характеризация выглядит неотличимо от честной.
            res["field_match"] = abs(res["B0_G"] - b0_guess) <= bracket_g
            if res["field_match"]:
                delta_note = ""
                if res.get("Delta_G") is None and res.get("Delta_G_from_gamma_bw") is not None:
                    delta_note = (f" (Delta_G не измерена напрямую -- "
                                   f"пересчитана из Gamma_BW/EPSUM_bg, "
                                   f"см. Delta_G_from_gamma_bw)")
                print(f"          -> B0={res['B0_G']:.6f} Гс, "
                      f"Delta={res.get('Delta_G')}, a_bg={res.get('a_bg_re')}"
                      f"{delta_note}")
                if res.get("Delta_G_from_gamma_bw") is not None:
                    print(f"          -> Delta_G_from_gamma_bw="
                          f"{res['Delta_G_from_gamma_bw']:.6g} Гс "
                          f"(Hutson NJP 9,152(2007) Eq.9, из "
                          f"Gamma_BW_G={res.get('Gamma_BW_G')})")
            else:
                print(f"          [WARNING] LOCPOL сошёлся на B0="
                      f"{res['B0_G']:.6f} Гс -- это ВНЕ исходного окна "
                      f"{b0_guess - bracket_g:.3f}-{b0_guess + bracket_g:.3f} Гс "
                      f"вокруг field.exe-оценки для state {row.state}. Похоже, "
                      f"экстраполяция 'убежала' на ДРУГОЙ резонанс (возможно, "
                      f"этот bound-уровень слабо связан с открытым каналом и "
                      f"не даёт наблюдаемого полюса рядом с {b0_guess:.6f} Гс) "
                      f"-- НЕ используйте Delta_G/a_bg этой строки как "
                      f"характеристику резонанса при {b0_guess:.6f} Гс, "
                      f"см. колонку field_match=False", file=sys.stderr)
        else:
            print(f"          -> НЕ охарактеризован (status={res.get('status')})")
        results.append(res)
    return results


def build_field_char_summary(results: list[dict], meta: RunMetadata) -> pd.DataFrame:
    """Как build_summary_table(), но для результатов
    characterize_field_resonances() -- добавляет колонки:
      state -- номер состояния из field.exe;
      field_G_field_driver -- B0 ДО уточнения через LOCPOL/IFCONV;
      field_match -- False означает, что LOCPOL сошёлся ВНЕ окна
        b0_guess +/- bracket_g (экстраполяция "убежала" на другой резонанс,
        см. предупреждение в characterize_field_resonances) -- ПРОВЕРЯЙТЕ
        эту колонку перед тем, как использовать Delta_G/a_bg строки: при
        field_match=False они относятся к СЛУЧАЙНО найденному резонансу,
        а не к тому, что предполагал field.exe для этого state (None -- для
        status != 'ok', сравнивать было не с чем)."""
    base = build_summary_table(results, meta)
    base.insert(0, "field_match", [r.get("field_match") for r in results])
    base.insert(0, "field_G_field_driver",
                [r.get("field_G_field_driver") for r in results])
    base.insert(0, "state", [r.get("state") for r in results])
    return base


def characterize_field_pipeline(cfg_field: Config, cfg_scatter: Config,
                                  fldmin: float, fldmax: float, monqn: str,
                                  bracket_g: float, label: str, meta: RunMetadata,
                                  csv_out: Path,
                                  l_filter: Optional[int] = 0) -> pd.DataFrame:
    """Полный мост field.exe -> LOCPOL/IFCONV: находит все B0 в [fldmin,fldmax]
    (Шаг 1, field.exe), затем характеризует КАЖДЫЙ (Шаг 2, molscat + IFCONV),
    получая Delta_G (=dB, ширина резонанса по полю) и a_bg. Пишет итоговый
    CSV и возвращает DataFrame."""
    output_dir = cfg_field.output_dir
    field_csv = output_dir / f"{label}_field_resonances.csv"

    print("=== Шаг 1/2: поиск точных B0 через field.exe ===")
    field_df = find_resonances_via_field(cfg_field, fldmin, fldmax, monqn,
                                          label, field_csv)

    if field_df.empty:
        print("\n[Шаг 2/2 пропущен] Характеризовать нечего -- field.exe не "
              "нашёл ни одного резонанса в этом окне.")
        summary = pd.DataFrame(columns=[
            "state", "field_G_field_driver", "field_match", "B0_G", "Delta_G",
            "a_bg_re", "a_bg_im", "a_res_re", "a_res_im", "Gamma_inel_G",
            "Gamma_BW_G", "EPSUM_bg_over_pi", "Delta_G_from_gamma_bw",
            "ifconv_used", "predicted_step_G",
            "channel", "L", "field_center_guess", "elapsed_sec",
            "n_molscat_calls", "status", "tensor_terms", "hyperfine_zeeman",
            "lmax", "dipdip_scale", "v2_scale", "run_label",
        ])
        summary.to_csv(csv_out, index=False)
        return summary

    print(f"\n=== Шаг 2/2: характеризация {len(field_df)} найденных "
          f"состояний (ICHAN + LOCPOL/IFCONV) ===")
    results = characterize_field_resonances(cfg_scatter, field_df, monqn,
                                             bracket_g, label, l_filter=l_filter)

    summary = build_field_char_summary(results, meta)
    summary.to_csv(csv_out, index=False)

    n_ok = (summary.status == "ok").sum()
    print(f"\n=== ГОТОВО: {n_ok}/{len(summary)} резонансов охарактеризовано, "
          f"записано в {csv_out} ===")
    return summary


# =========================================================================
# 7. ВЕРХНЕУРОВНЕВЫЙ ПАЙПЛАЙН
# =========================================================================

def run_full_pipeline(cfg: Config, meta: RunMetadata, fmin: float, fmax: float,
                       dfield: float, csv_out: Path, l_filter: Optional[int] = 0,
                       min_separation_G: float = 15.0, window_mult: float = 0.25) -> pd.DataFrame:
    t_start = time.perf_counter()
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    coarse_dir = cfg.output_dir / f"{meta.label}_coarse"
    coarse_dir.mkdir(parents=True, exist_ok=True)

    print(f"[1/6] Генерация input грубого скана {fmin}-{fmax} Гс, шаг {dfield} Гс")
    coarse_input = coarse_dir / "coarse.input"
    generate_coarse_scan_input(
        cfg, fmin=fmin, fmax=fmax, dfield=dfield,
        label=f"{meta.label} coarse scan", out_path=coarse_input,
    )

    print("[2/6] Запуск molscat (грубый скан)...")
    coarse_out, coarse_elapsed = run_molscat(cfg, coarse_input)
    print(f"       -> заняло {_fmt_duration(coarse_elapsed)}")

    print("[2/6] Парсинг вывода (ВСЕ каналы, без фильтра по ICHAN)...")
    df = parse_coarse_scan(coarse_out, target_ichan=None)
    print(f"       -> {len(df)} строк (полей x каналов), "
          f"{df.channel.nunique() if not df.empty else 0} уникальных каналов")

    print(f"[3/6] Автодетекция кандидатов (L={l_filter}, "
          f"min_separation_G={min_separation_G}, window_mult={window_mult})...")
    candidates = detect_candidates(df, l_filter=l_filter,
                                    min_separation_G=min_separation_G,
                                    window_mult=window_mult)
    print(f"       -> найдено {len(candidates)} кандидатов: "
          f"{[(round(c['field_center'],1), 'chan='+str(c['channel'])) for c in candidates]}")

    if not candidates:
        print("       -> кандидатов не найдено -- характеризация и сборка "
              "таблицы пропущены (это НЕ обязательно ошибка: возможно, в "
              "этом диапазоне поля действительно нет резонансов -- либо "
              "--min-separation-g/--window-mult/--l-filter настроены не под "
              "этот скан, см. docstring detect_candidates(); проверьте "
              f"{coarse_out} глазами, прежде чем сужать диапазон дальше)")
        summary = pd.DataFrame(columns=[
            "B0_G", "Delta_G", "a_bg_re", "a_bg_im", "a_res_re", "a_res_im",
            "Gamma_inel_G", "Gamma_BW_G", "EPSUM_bg_over_pi", "ifconv_used",
            "predicted_step_G", "channel", "L", "field_center_guess",
            "elapsed_sec", "n_molscat_calls", "status", "tensor_terms",
            "hyperfine_zeeman", "lmax", "dipdip_scale", "v2_scale",
            "run_label",
        ])
        summary.to_csv(csv_out, index=False)
        print(f"       -> записан пустой CSV (0 кандидатов) в {csv_out}")
        total_elapsed = time.perf_counter() - t_start
        print(f"\n=== ГОТОВО за {_fmt_duration(total_elapsed)} "
              f"(скан: {_fmt_duration(coarse_elapsed)}, кандидатов не "
              f"найдено, характеризация не запускалась) ===")
        return summary

    print("[4-5/6] Характеризация каждого кандидата (IFCONV с эскалацией)...")
    results = []
    cand_durations = []
    n_cand = len(candidates)
    for i, c in enumerate(candidates):
        t_cand = time.perf_counter()
        print(f"       -> кандидат {i+1}/{n_cand} @ {c['field_center']:.1f} Гс, "
              f"канал {c['channel']} ({c['reason']})")
        run_id = f"cand{i+1}_{c['field_center']:.0f}G"
        res = characterize_candidate(cfg, c, run_id)
        dt = time.perf_counter() - t_cand
        cand_durations.append(dt)
        avg_so_far = sum(cand_durations) / len(cand_durations)
        eta = avg_so_far * (n_cand - i - 1)
        print(f"          статус: {res.get('status')}  "
              f"[{_fmt_duration(dt)}, {res.get('n_molscat_calls', '?')} запусков molscat]  "
              f"осталось ~{_fmt_duration(eta)} ({n_cand-i-1} кандидатов)")
        results.append(res)

    print("[6/6] Сборка итоговой таблицы...")
    summary = build_summary_table(results, meta)
    n_before = (summary.status == "ok").sum()
    summary = dedupe_by_B0(summary)
    n_after = (summary.status == "ok").sum()
    if n_after < n_before:
        print(f"       -> убрано {n_before - n_after} дублей (один и тот же "
              f"резонанс, найденный с двух соседних узлов сетки)")
    summary.to_csv(csv_out, index=False)
    print(f"       -> записано в {csv_out}")

    total_elapsed = time.perf_counter() - t_start
    n_ok = (summary.status == "ok").sum()
    n_failed = (summary.status != "ok").sum()
    print(f"\n=== ГОТОВО за {_fmt_duration(total_elapsed)} "
          f"(скан: {_fmt_duration(coarse_elapsed)}, "
          f"характеризация {n_cand} кандидатов: "
          f"{_fmt_duration(sum(cand_durations))}) ===")
    print(f"    {n_ok} сошлось, {n_failed} нет. "
          f"В среднем {_fmt_duration(sum(cand_durations)/n_cand) if n_cand else '0s'} на кандидата.")
    print(f"    Все .input/.out лежат в {cfg.output_dir}")
    return summary


# =========================================================================
# CLI
# =========================================================================

def _selftest_parse(out_file: str, ifconv: int = 1):
    """Диагностика: прогнать парсер на реальном .out и показать, что нашлось."""
    p = Path(out_file)
    df = parse_coarse_scan(p)
    print(f"Найдено {len(df)} строк таблицы K-DEPENDENT SCATTERING LENGTHS.")
    if not df.empty:
        print(df.head(20).to_string())
    else:
        print("НИЧЕГО не найдено (это нормально, если файл -- вывод IFCONV-"
              "характеризации, а не грубого скана; таблица K-DEPENDENT там "
              "тоже печатается на каждой итерации, так что пустой результат "
              "по-настоящему подозрителен только для coarse-scan файлов).")

    result = parse_ifconv_block(p, ifconv)
    print(f"\nparse_ifconv_block(ifconv={ifconv}) -> {result}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_full = sub.add_parser("run", help="полный конвейер: coarse scan -> кандидаты -> характеризация -> CSV")
    p_full.add_argument("--exe", required=True, type=Path, help="путь к molscat-Tm2.exe")
    p_full.add_argument("--work-dir", required=True, type=Path,
                         help="Рабочая директория для запуска molscat-Tm2.exe "
                              "(там, где бинарник и его зависимости, обычно "
                              "source_code). НЕ используется для хранения "
                              ".input/.out -- см. --output-dir.")
    p_full.add_argument("--output-dir", type=Path, default=None,
                         help="Куда складывать все сгенерированные .input/.out "
                              "(отдельно от --work-dir, чтобы не засорять "
                              "source_code). По умолчанию: <work-dir>/molscat_runs. "
                              "Внутри создаётся подпапка на каждый кандидат.")
    p_full.add_argument("--fmin", type=float, default=0.0)
    p_full.add_argument("--fmax", type=float, default=2000.0)
    p_full.add_argument("--dfield", type=float, default=20.0)
    p_full.add_argument("--lmax", type=int, default=4)
    p_full.add_argument("--jtot", type=int, default=-12,
                         help="JTOT = 2*MTOT (см. базис в base9-Tm_Tm_AI.f). "
                              "Дефолт -12 (обычный входной канал mF=-6 на "
                              "пару). Для полностью растянутого канала "
                              "mF=-8,-8 используйте -16 (как в MONQN "
                              "'8,-8,8,-8').")
    p_full.add_argument("--ibfix", type=int, default=2,
                         help="Симметричный блок (чётность). Дефолт 2.")
    p_full.add_argument("--jstep", type=int, default=2,
                         help="Дефолт 2 (см. Config.jstep).")
    p_full.add_argument("--label", default="V0V2_run1")
    p_full.add_argument("--dipdip-scale", type=float, default=1.0,
                         dest="dipdip_scale",
                         help="Множитель перед D^(2)_2 (магнитный "
                              "диполь-диполь, k=2,i=2, [j1⊗j2]_2q -- см. "
                              "DIPBLK/TENSX в base9-Tm_Tm_AI.f). 1.0 (дефолт) "
                              "= реальный физический коэффициент, 0.0 = "
                              "член выключен (V^(2)_2(R)=0), любое другое "
                              "значение -- искусственное масштабирование "
                              "для сравнительных прогонов. Требует бинарник, "
                              "собранный с MXLAM=3 в base9-Tm_Tm_AI.f "
                              "(POTIN9) -- шаблон &POTL теперь ВСЕГДА пишет "
                              "3 блока; со старым бинарником под MXLAM=2 "
                              "molscat завершится с ошибкой чтения &POTL.")
    p_full.add_argument("--v2-scale", type=float, default=1.0,
                         dest="v2_scale",
                         help="Множитель перед коэффициентом связи блока "
                              "k=2,i=1 (анизотропный V^(1)_2(R), Table 5 -- "
                              "см. V2SCALE в &BASIS9 / base9-Tm_Tm_AI.f). "
                              "Физически эквивалентно масштабированию "
                              "V^(1)_2(R) целиком для всех R. 1.0 (дефолт) "
                              "= реальная физика, 0.0 = анизотропия k=2,i=1 "
                              "выключена, любое другое значение -- "
                              "искусственное масштабирование для "
                              "сравнительных прогонов. Требует бинарник, "
                              "собранный с V2SCALE в NAMELIST /BASIS9/ "
                              "base9-Tm_Tm_AI.f -- со старым бинарником "
                              "molscat завершится с ошибкой чтения &BASIS9.")
    p_full.add_argument("--l-filter", type=int, default=0,
                         help="Партиальная волна L, по которой ищем кандидатов "
                              "в грубом скане (0 = s-волна, физически стандартная "
                              "длина рассеяния; -1 = отключить фильтр и искать по "
                              "всем L -- НЕ рекомендуется, т.к. L>0 даёт другую "
                              "размерность величины и ложные 'выбросы'). "
                              "ICHAN для характеризации каждого кандидата теперь "
                              "подбирается автоматически, вручную задавать не нужно.")
    p_full.add_argument("--min-separation-g", type=float, default=15.0,
                         dest="min_separation_g",
                         help="Минимальное расстояние (Гс) между кандидатами, "
                              "при котором они схлопываются в одного (иначе "
                              "считаются разными резонансами). Дефолт 15.0 "
                              "откалиброван под грубый шаг 20 Гс; при мелком "
                              "шаге (напр. пилотный скан на узком окне с "
                              "dfield=0.2) уменьшайте пропорционально -- иначе "
                              "весь диапазон может схлопнуться в 1 кандидата.")
    p_full.add_argument("--window-mult", type=float, default=0.25,
                         help="Множитель окна характеризации кандидата: "
                              "field_lo/hi = field_center +/- dfield*window_mult. "
                              "Дефолт 0.25 (т.е. +/-dfield/4) откалиброван под "
                              "шаг 20 Гс (окно ~5 Гс). При мелком dfield "
                              "увеличивайте (напр. window_mult~5-10 при "
                              "dfield=0.2, чтобы окно оставалось ~1-2 Гс -- "
                              "заметно шире ожидаемой Delta, но не настолько "
                              "широким, чтобы зацепить соседний резонанс).")
    p_full.add_argument("--csv-out", type=Path, default=None,
                         help="По умолчанию: <output-dir>/<label>_summary.csv")
    p_full.add_argument("--tensor-terms", default="k0i1,k2i1,k2i2",
                         help="Какие (k,i)-члены V(R) реально включены в этот "
                              "прогон (см. Table 1, Tiesinga et al. NJP 23, "
                              "085007) -- ТОЛЬКО метаданные для CSV, поставьте "
                              "в соответствие --dipdip-scale вручную (напр. "
                              "'k0i1,k2i1' при --dipdip-scale 0, "
                              "'k0i1,k2i1,k2i2' при --dipdip-scale 1) -- этот "
                              "флаг физику не переключает, физику переключает "
                              "--dipdip-scale.")
    p_full.add_argument("--monqn", default=None,
                         help="Задать входной канал через физические "
                              "квантовые числа F,mF обоих атомов вместо "
                              "IREF=1 (использует THRSH9/Брейт-Раби в "
                              "base9-Tm_Tm_AI.f -- решает проблему "
                              "вырождения EREF при B->0 и необходимость "
                              "переподбирать ICHAN для каждого кандидата "
                              "по отдельности). Формат: 4 ДОУБЛЕННЫХ числа "
                              "'2F1,2mF1,2F2,2mF2' через запятую, напр. "
                              "'8,-8,8,-8' для F1=F2=4,mF1=mF2=-4 (обычный "
                              "полностью растянутый входной канал). По "
                              "умолчанию не задано -- старое поведение "
                              "(IREF=1). ICHAN всё равно нужно передавать "
                              "отдельно (через --l-filter/candidate "
                              "channel), MONQN влияет только на EREF.")
    p_full.add_argument("--no-hyperfine-zeeman", action="store_true",
                         help="Пометить в CSV, что одноатомные (HFSPLA/GSA) "
                              "члены выключены -- переключатель для этого в "
                              "base9-Tm_Tm_AI.f пока не реализован (Этап 4).")

    p_test = sub.add_parser("selftest-parse", help="диагностика регэкспов на реальном .out")
    p_test.add_argument("out_file")
    p_test.add_argument("ifconv", type=int, nargs="?", default=1,
                         help="IFCONV, с которым был запущен этот .out (1/2/3); "
                              "влияет на то, какой формат разбора (130 vs 106) "
                              "использовать. По умолчанию 1.")

    p_char = sub.add_parser("characterize",
                             help="точечная характеризация ОДНОГО кандидата "
                                  "(без полного скана) -- удобно для перезапуска "
                                  "после фикса/для проблемных резонансов")
    p_char.add_argument("--exe", required=True, type=Path)
    p_char.add_argument("--work-dir", required=True, type=Path)
    p_char.add_argument("--output-dir", type=Path, default=None)
    p_char.add_argument("--fmin", type=float, required=True,
                         help="Нижняя граница бракетирующего окна, Гс")
    p_char.add_argument("--fmax", type=float, required=True,
                         help="Верхняя граница бракетирующего окна, Гс")
    p_char.add_argument("--ichan", type=int, required=True,
                         help="ICHAN -- фиксированный номер входного канала "
                              "в базисе (см. таблицу K-DEPENDENT... грубого "
                              "скана для этого поля)")
    p_char.add_argument("--l", type=int, default=0, dest="l_value",
                         help="L входного канала (только для справки в CSV)")
    p_char.add_argument("--lmax", type=int, default=4)
    p_char.add_argument("--jtot", type=int, default=-12,
                         help="JTOT = 2*MTOT (см. базис в base9-Tm_Tm_AI.f). "
                              "Дефолт -12. Для полностью растянутого канала "
                              "mF=-8,-8 используйте -16.")
    p_char.add_argument("--ibfix", type=int, default=2,
                         help="Симметричный блок (чётность). Дефолт 2.")
    p_char.add_argument("--jstep", type=int, default=2,
                         help="Дефолт 2 (см. Config.jstep).")
    p_char.add_argument("--label", default="manual_cand")
    p_char.add_argument("--monqn", default=None,
                         help="См. help в 'run' -- '2F1,2mF1,2F2,2mF2', "
                              "напр. '8,-8,8,-8'.")
    p_char.add_argument("--dipdip-scale", type=float, default=1.0,
                         dest="dipdip_scale",
                         help="См. help в 'run' -- множитель перед D^(2)_2.")
    p_char.add_argument("--v2-scale", type=float, default=1.0,
                         dest="v2_scale",
                         help="См. help в 'run' -- множитель перед V^(1)_2(R).")

    p_find = sub.add_parser("find-channel",
                             help="найти номер ICHAN, соответствующий "
                                  "заданным (F,mF) обоих атомов, на "
                                  "конкретном поле (см. THRSH9 в "
                                  "base9-Tm_Tm_AI.f -- сам он номер "
                                  "канала не знает, только EREF)")
    p_find.add_argument("--exe", required=True, type=Path)
    p_find.add_argument("--work-dir", required=True, type=Path)
    p_find.add_argument("--output-dir", type=Path, default=None)
    p_find.add_argument("--monqn", required=True,
                         help="'2F1,2mF1,2F2,2mF2', напр. '8,-8,8,-8'.")
    p_find.add_argument("--field", type=float, required=True,
                         help="Поле (Гс), на котором ищем номер канала. "
                              "ИЗБЕГАЙТЕ B=0 -- разные (F,mF) могут быть "
                              "там точно вырождены (см. обсуждение "
                              "вырождения при B->0); возьмите любое "
                              "малое ненулевое поле, например 0.5-1 Гс, "
                              "или сразу поле интересующего резонанса -- "
                              "номер канала не должен меняться внутри "
                              "одного JTOT/IBFIX сектора (адиабатическое "
                              "слежение), кроме как в точках истинного "
                              "пересечения уровней.")
    p_find.add_argument("--jtot", type=int, default=-12)
    p_find.add_argument("--ibfix", type=int, default=2)
    p_find.add_argument("--jstep", type=int, default=2)
    p_find.add_argument("--lmax", type=int, default=4)
    p_find.add_argument("--label", default="find_channel")
    p_find.add_argument("--dipdip-scale", type=float, default=1.0,
                         dest="dipdip_scale",
                         help="См. help в 'run' -- множитель перед D^(2)_2.")
    p_find.add_argument("--v2-scale", type=float, default=1.0,
                         dest="v2_scale",
                         help="См. help в 'run' -- множитель перед V^(1)_2(R).")
    p_find.add_argument("--tol-cm1", type=float, default=1.0e-8, dest="tol_cm1",
                         help="Допуск (см⁻¹) при сравнении EREF с "
                              "энергиями из таблицы THRESHOLDS.")
    p_find.add_argument("--l-filter", type=int, default=0, dest="l_filter",
                         help="Партциальная волна L для разрешения "
                              "неоднозначности ICHAN (см. find_channel_for_monqn: "
                              "при LMAX>=2 в таблице THRESHOLDS ВСЕГДА несколько "
                              "строк с одинаковой энергией порога -- по одной на "
                              "каждый L, т.к. энергия порога от L не зависит; "
                              "это не 'настоящее' вырождение уровней). Дефолт "
                              "0 (s-волна, стандартный выбор). -1 отключает "
                              "фильтр (вернуться к старому поведению -- "
                              "'ambiguous' при LMAX>=2 почти всегда).")

    p_field = sub.add_parser("find-resonances",
                              help="найти позиции резонансов НАПРЯМУЮ через "
                                   "field.exe (без грубого скана molscat) -- "
                                   "ищет поле(я) в [--fmin,--fmax], где "
                                   "связанный уровень пересекает порог. "
                                   "ЛУЧШЕ подходит для ситуации с одним "
                                   "открытым каналом, чем 'run'.")
    p_field.add_argument("--exe", required=True, type=Path,
                          help="путь к field-Tm2.exe (НЕ molscat-Tm2.exe!)")
    p_field.add_argument("--work-dir", required=True, type=Path)
    p_field.add_argument("--output-dir", type=Path, default=None)
    p_field.add_argument("--fmin", type=float, required=True, dest="fldmin",
                          help="Нижняя граница окна поиска, Гс")
    p_field.add_argument("--fmax", type=float, required=True, dest="fldmax",
                          help="Верхняя граница окна поиска, Гс")
    p_field.add_argument("--monqn", required=True,
                          help="ОБЯЗАТЕЛЕН (в отличие от 'run'/'characterize') "
                               "-- '2F1,2mF1,2F2,2mF2', напр. '8,-8,8,-8' "
                               "(F1=F2=4,mF1=mF2=-4, полностью растянутый "
                               "входной канал для Tm2, I=1/2 -- см. "
                               "ABS(2F-ISA)==1 в THRSH9 base9-Tm_Tm_AI.f). "
                               "Без него не определён порог, относительно "
                               "которого ищется пересечение.")
    p_field.add_argument("--jtot", type=int, default=-12)
    p_field.add_argument("--ibfix", type=int, default=2)
    p_field.add_argument("--jstep", type=int, default=2)
    p_field.add_argument("--lmax", type=int, default=4)
    p_field.add_argument("--ured", type=float, default=84.467109,
                          help="Приведённая масса пары (а.е.м.). Дефолт -- "
                               "физическая масса пары 169Tm.")
    p_field.add_argument("--label", default="field_search")
    p_field.add_argument("--dipdip-scale", type=float, default=1.0,
                          dest="dipdip_scale",
                          help="См. help в 'run' -- множитель перед D^(2)_2.")
    p_field.add_argument("--v2-scale", type=float, default=1.0,
                          dest="v2_scale",
                          help="См. help в 'run' -- множитель перед V^(1)_2(R).")
    p_field.add_argument("--csv-out", type=Path, default=None,
                          help="По умолчанию: <output-dir>/<label>_resonances.csv")

    p_charfield = sub.add_parser(
        "characterize-field",
        help="полная эластичная процедура Frye&Hutson (PRA 96, 042705, Sec. II): "
             "field.exe находит точные B0, затем для каждого автоматически "
             "определяется ICHAN и запускается LOCPOL/IFCONV -- даёт B0, "
             "Delta_G (=dB, ширина резонанса по полю) и a_bg одним прогоном "
             "на резонанс, без ручного подбора окон/каналов.")
    p_charfield.add_argument("--field-exe", required=True, type=Path,
                              dest="field_exe",
                              help="путь к field-Tm2.exe (Шаг 1: поиск B0)")
    p_charfield.add_argument("--molscat-exe", required=True, type=Path,
                              dest="molscat_exe",
                              help="путь к molscat-Tm2.exe (Шаг 2: "
                                   "характеризация через IFCONV/LOCPOL)")
    p_charfield.add_argument("--work-dir", required=True, type=Path)
    p_charfield.add_argument("--output-dir", type=Path, default=None)
    p_charfield.add_argument("--fmin", type=float, required=True, dest="fldmin",
                              help="Нижняя граница окна поиска B0, Гс")
    p_charfield.add_argument("--fmax", type=float, required=True, dest="fldmax",
                              help="Верхняя граница окна поиска B0, Гс")
    p_charfield.add_argument("--monqn", required=True,
                              help="'2F1,2mF1,2F2,2mF2', напр. '8,-8,8,-8' "
                                   "(F1=F2=4,mF1=mF2=-4, полностью растянутый "
                                   "входной канал для Tm2, I=1/2) -- "
                                   "используется И для field.exe, И для "
                                   "последующего определения ICHAN/IFCONV "
                                   "(должен быть ОДИН и тот же канал на "
                                   "обоих шагах).")
    p_charfield.add_argument("--bracket-g", type=float, default=0.5,
                              dest="bracket_g",
                              help="Начальная полуширина бракетирующего окна "
                                   "(Гс) вокруг каждого B0 от field.exe для "
                                   "старта LOCPOL/IFCONV -- см. Sec. II "
                                   "статьи (δB=0.2 Гс между начальными "
                                   "точками). B0 здесь уже точен (в отличие "
                                   "от кандидатов грубого скана), так что "
                                   "окно нужно только чтобы дать LOCPOL "
                                   "стартовые точки; если не подойдёт -- "
                                   "characterize_candidate сам его сузит. "
                                   "Дефолт 0.5 Гс; уменьшайте, если "
                                   "резонансы расположены ближе друг к "
                                   "другу, чем 1 Гс.")
    p_charfield.add_argument("--jtot", type=int, default=-12)
    p_charfield.add_argument("--ibfix", type=int, default=2)
    p_charfield.add_argument("--jstep", type=int, default=2)
    p_charfield.add_argument("--lmax", type=int, default=4)
    p_charfield.add_argument("--ured", type=float, default=84.467109,
                              help="Приведённая масса пары (а.е.м.). Дефолт "
                                   "-- физическая масса пары 169Tm.")
    p_charfield.add_argument("--label", default="charfield")
    p_charfield.add_argument("--dipdip-scale", type=float, default=1.0,
                              dest="dipdip_scale",
                              help="См. help в 'run' -- множитель перед D^(2)_2.")
    p_charfield.add_argument("--v2-scale", type=float, default=1.0,
                              dest="v2_scale",
                              help="См. help в 'run' -- множитель перед V^(1)_2(R).")
    p_charfield.add_argument("--tensor-terms", default="k0i1,k2i1,k2i2",
                              help="См. help в 'run' -- только метаданные CSV.")
    p_charfield.add_argument("--no-hyperfine-zeeman", action="store_true",
                              help="См. help в 'run' -- только метаданные CSV.")
    p_charfield.add_argument("--l-filter", type=int, default=0, dest="l_filter",
                              help="См. find-channel --l-filter -- разрешает "
                                   "неизбежную L-неоднозначность ICHAN при "
                                   "LMAX>=2. Дефолт 0 (s-волна). -1 отключает "
                                   "фильтр.")
    p_charfield.add_argument("--csv-out", type=Path, default=None,
                              help="По умолчанию: <output-dir>/<label>_char_summary.csv")

    p_count = sub.add_parser("count-resonances",
                              help="БЫСТРО посчитать ЧИСЛО резонансов в окне "
                                   "[--fmin,--fmax] через разницу счёта узлов "
                                   "на FLDMIN/FLDMAX (то же, что field.exe сам "
                                   "печатает как 'SEEK N STATES', но с "
                                   "MXCALC=3 -- без дорогого поиска точных "
                                   "позиций, см. find-resonances для этого).")
    p_count.add_argument("--exe", required=True, type=Path,
                          help="путь к field-Tm2.exe (НЕ molscat-Tm2.exe!)")
    p_count.add_argument("--work-dir", required=True, type=Path)
    p_count.add_argument("--output-dir", type=Path, default=None)
    p_count.add_argument("--fmin", type=float, required=True, dest="fldmin",
                          help="Нижняя граница окна, Гс")
    p_count.add_argument("--fmax", type=float, required=True, dest="fldmax",
                          help="Верхняя граница окна, Гс")
    p_count.add_argument("--monqn", required=True,
                          help="ОБЯЗАТЕЛЕН -- '2F1,2mF1,2F2,2mF2', напр. "
                               "'8,-8,8,-8'. См. help в 'find-resonances'.")
    p_count.add_argument("--jtot", type=int, default=-12)
    p_count.add_argument("--ibfix", type=int, default=2)
    p_count.add_argument("--jstep", type=int, default=2)
    p_count.add_argument("--lmax", type=int, default=4)
    p_count.add_argument("--ured", type=float, default=84.467109,
                          help="Приведённая масса пары (а.е.м.). Дефолт -- "
                               "физическая масса пары 169Tm.")
    p_count.add_argument("--label", default="node_count")
    p_count.add_argument("--dipdip-scale", type=float, default=1.0,
                          dest="dipdip_scale",
                          help="См. help в 'run' -- множитель перед D^(2)_2.")
    p_count.add_argument("--v2-scale", type=float, default=1.0,
                          dest="v2_scale",
                          help="См. help в 'run' -- множитель перед V^(1)_2(R).")

    args = ap.parse_args()

    if args.cmd == "selftest-parse":
        _selftest_parse(args.out_file, args.ifconv)
        return

    if args.cmd == "characterize":
        args.work_dir.mkdir(parents=True, exist_ok=True)
        output_dir = args.output_dir or (args.work_dir / "molscat_runs")
        output_dir.mkdir(parents=True, exist_ok=True)
        cfg = Config(
            molscat_exe=args.exe, work_dir=args.work_dir, output_dir=output_dir,
            coarse_template=Path(""), lmax=args.lmax, monqn=args.monqn,
            jtot=args.jtot, ibfix=args.ibfix, jstep=args.jstep,
            dipdip_scale=args.dipdip_scale, v2_scale=args.v2_scale,
        )
        candidate = {
            "field_lo": args.fmin, "field_hi": args.fmax,
            "field_center": (args.fmin + args.fmax) / 2,
            "channel": args.ichan, "L": args.l_value,
            "reason": "ручной перезапуск (characterize)",
        }
        t0 = time.perf_counter()
        res = characterize_candidate(cfg, candidate, args.label)
        dt = time.perf_counter() - t0
        print(f"\n=== ГОТОВО за {_fmt_duration(dt)}, "
              f"{res.get('n_molscat_calls', '?')} запусков molscat ===")
        print(f"статус: {res.get('status')}")
        if res.get("status") == "ok":
            for k in ("B0_G", "Delta_G", "a_bg_re", "a_bg_im", "a_res_re",
                      "a_res_im", "Gamma_inel_G", "Gamma_BW_G",
                      "EPSUM_bg_over_pi", "ifconv_used"):
                if k in res:
                    print(f"  {k} = {res[k]}")
        else:
            print(f"  last_status: {res.get('last_status')}")
            if res.get("tentative"):
                print(f"  tentative: {res['tentative']}")
        print(f"  файлы: {cfg.output_dir / args.label}")
        return

    if args.cmd == "find-channel":
        args.work_dir.mkdir(parents=True, exist_ok=True)
        output_dir = args.output_dir or (args.work_dir / "molscat_runs")
        output_dir.mkdir(parents=True, exist_ok=True)
        cfg = Config(
            molscat_exe=args.exe, work_dir=args.work_dir, output_dir=output_dir,
            coarse_template=Path(""), lmax=args.lmax, jtot=args.jtot,
            ibfix=args.ibfix, jstep=args.jstep, dipdip_scale=args.dipdip_scale,
            v2_scale=args.v2_scale,
        )
        l_filter = None if args.l_filter < 0 else args.l_filter
        res = find_channel_for_monqn(cfg, args.monqn, args.field, args.label,
                                      tol_cm1=args.tol_cm1, l_filter=l_filter)
        print(f"\nMONQN = {args.monqn}, поле = {args.field} Гс")
        if res["status"] == "ok":
            print(f"=== НАЙДЕН КАНАЛ: ICHAN = {res['channel']} "
                  f"(L = {res['L']}) ===")
            print(f"  energy_cm1  = {res['energy_cm1']!r}")
            print(f"  EREF (THRSH9) = {res['reference_energy_cm1']!r}")
            print(f"  diff        = {res['diff_cm1']:.3e} см-1")
        elif res["status"] == "ambiguous":
            print(f"=== НЕОДНОЗНАЧНО: {res['reason']} ===")
            for c in res["candidates"]:
                print(f"  ICHAN={c['index']}  L={c['L']}  "
                      f"energy_cm1={c['energy_cm1']!r}")
        else:
            print(f"=== НЕ НАЙДЕНО: {res.get('reason')} ===")
            if "all_thresholds" in res:
                print(f"  EREF (THRSH9) = {res['reference_energy_cm1']!r}")
                print("  все пороги в таблице:")
                for t in res["all_thresholds"]:
                    print(f"    ICHAN={t['index']}  L={t['L']}  "
                          f"energy_cm1={t['energy_cm1']!r}")
        print(f"  файл: {res.get('out_file')}")
        return

    if args.cmd == "find-resonances":
        args.work_dir.mkdir(parents=True, exist_ok=True)
        output_dir = args.output_dir or (args.work_dir / "molscat_runs")
        output_dir.mkdir(parents=True, exist_ok=True)
        csv_out = args.csv_out or (output_dir / f"{args.label}_resonances.csv")
        cfg = Config(
            molscat_exe=args.exe, work_dir=args.work_dir, output_dir=output_dir,
            coarse_template=Path(""), lmax=args.lmax, jtot=args.jtot,
            ibfix=args.ibfix, jstep=args.jstep, dipdip_scale=args.dipdip_scale,
            v2_scale=args.v2_scale, ured=args.ured,
        )
        find_resonances_via_field(cfg, args.fldmin, args.fldmax, args.monqn,
                                   args.label, csv_out)
        return

    if args.cmd == "characterize-field":
        args.work_dir.mkdir(parents=True, exist_ok=True)
        output_dir = args.output_dir or (args.work_dir / "molscat_runs")
        output_dir.mkdir(parents=True, exist_ok=True)
        csv_out = args.csv_out or (output_dir / f"{args.label}_char_summary.csv")

        cfg_field = Config(
            molscat_exe=args.field_exe, work_dir=args.work_dir,
            output_dir=output_dir, coarse_template=Path(""), lmax=args.lmax,
            jtot=args.jtot, ibfix=args.ibfix, jstep=args.jstep,
            dipdip_scale=args.dipdip_scale, v2_scale=args.v2_scale,
            ured=args.ured,
        )
        cfg_scatter = dataclasses.replace(
            cfg_field, molscat_exe=args.molscat_exe, monqn=args.monqn,
        )
        meta = RunMetadata(
            tensor_terms=args.tensor_terms,
            hyperfine_zeeman=not args.no_hyperfine_zeeman,
            lmax=args.lmax, label=args.label,
            dipdip_scale=args.dipdip_scale, v2_scale=args.v2_scale,
        )
        l_filter = None if args.l_filter < 0 else args.l_filter
        characterize_field_pipeline(
            cfg_field, cfg_scatter, args.fldmin, args.fldmax, args.monqn,
            args.bracket_g, args.label, meta, csv_out, l_filter=l_filter,
        )
        return

    if args.cmd == "count-resonances":
        args.work_dir.mkdir(parents=True, exist_ok=True)
        output_dir = args.output_dir or (args.work_dir / "molscat_runs")
        output_dir.mkdir(parents=True, exist_ok=True)
        cfg = Config(
            molscat_exe=args.exe, work_dir=args.work_dir, output_dir=output_dir,
            coarse_template=Path(""), lmax=args.lmax, jtot=args.jtot,
            ibfix=args.ibfix, jstep=args.jstep, dipdip_scale=args.dipdip_scale,
            v2_scale=args.v2_scale, ured=args.ured,
        )
        res = count_resonances_via_field(cfg, args.fldmin, args.fldmax,
                                          args.monqn, args.label)
        print(f"\n=== N резонансов в [{args.fldmin},{args.fldmax}] Гс: "
              f"{res.get('n_resonances')} (status={res.get('status')}) ===")
        return

    if args.cmd == "run":
        args.work_dir.mkdir(parents=True, exist_ok=True)
        output_dir = args.output_dir or (args.work_dir / "molscat_runs")
        output_dir.mkdir(parents=True, exist_ok=True)
        csv_out = args.csv_out or (output_dir / f"{args.label}_summary.csv")

        cfg = Config(
            molscat_exe=args.exe, work_dir=args.work_dir, output_dir=output_dir,
            coarse_template=Path(""),  # не используется, шаблон встроен
            lmax=args.lmax, monqn=args.monqn,
            jtot=args.jtot, ibfix=args.ibfix, jstep=args.jstep,
            dipdip_scale=args.dipdip_scale, v2_scale=args.v2_scale,
        )
        meta = RunMetadata(
            tensor_terms=args.tensor_terms,
            hyperfine_zeeman=not args.no_hyperfine_zeeman,
            lmax=args.lmax, label=args.label,
            dipdip_scale=args.dipdip_scale, v2_scale=args.v2_scale,
        )
        l_filter = None if args.l_filter < 0 else args.l_filter
        run_full_pipeline(cfg, meta, args.fmin, args.fmax, args.dfield,
                           csv_out, l_filter=l_filter,
                           min_separation_G=args.min_separation_g,
                           window_mult=args.window_mult)


if __name__ == "__main__":
    main()
