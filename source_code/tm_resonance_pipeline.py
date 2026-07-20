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
    В base9-Tm_Tm_AI.f реализованы и верифицированы ТОЛЬКО (0,1) и (2,1)."""
    tensor_terms: str = "k0i1,k2i1"  # список включённых (k,i)-членов V(R),
                                       # запятая между ними; текущий базовый
                                       # прогон = только то, что верифицировано
    hyperfine_zeeman: bool = True     # одноатомные члены (HFSPLA/GSA в
                                       # BAS9IN) -- НЕ часть разложения V(R),
                                       # отдельная ось; сейчас всегда включены
                                       # (переключатель "выкл" -- Этап 4, ещё
                                       # не реализован)
    lmax: int = 4                     # LMAX БАЗИСА (число партиальных волн
                                       # в разложении по l), это НЕ тензорный
                                       # ранг k -- не путать при чтении CSV!
    label: str = "k0i1_k2i1_only"


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
    LMAX   =  {lmax},
 /

 &POTL
    MXLAM  = 2, LAMBDA =  {lambda_terms},
                NTERM  = -1, -1,
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
    LMAX   =  {lmax},
 /

 &POTL
    MXLAM  = 2, LAMBDA =  {lambda_terms},
                NTERM  = -1, -1,
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


def generate_coarse_scan_input(cfg: Config, fmin: float, fmax: float,
                                 dfield: float, label: str,
                                 out_path: Path) -> Path:
    iref, monqn_clause = _iref_and_monqn_clause(cfg.monqn)
    text = COARSE_TEMPLATE.format(
        label=label, ured=cfg.ured, jtot=cfg.jtot, ibfix=cfg.ibfix,
        jstep=cfg.jstep, energy=_fortran_e(cfg.energy_K), fmin=fmin, fmax=fmax,
        dfield=dfield, ichan=cfg.ichan_guess, lmax=cfg.lmax,
        lambda_terms=cfg.lambda_terms, iref=iref, monqn_clause=monqn_clause,
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
        iref=iref, monqn_clause=monqn_clause,
    )
    out_path.write_text(text)
    return out_path


# =========================================================================
# 2. ЗАПУСК MOLSCAT
# =========================================================================

def run_molscat(cfg: Config, input_path: Path, timeout: int = 1800) -> tuple[Path, float]:
    """Запускает molscat-Tm2.exe < input_path > input_path.with_suffix('.out').
    input_path/output .out живут в cfg.output_dir; сам процесс запускается с
    cwd=cfg.work_dir (там, где бинарник и его зависимости -- НЕ то же самое,
    что output_dir, см. Config).
    Возвращает (путь к .out, время выполнения в секундах). Не бросает
    исключение при ненулевом коде возврата — molscat иногда возвращает его
    даже при штатных сообщениях типа NOPEN CHANGED; решение о фатальности
    принимает вызывающий код по содержимому .out."""
    out_path = input_path.with_suffix(".out")
    t0 = time.perf_counter()
    with open(input_path, "r") as fin, open(out_path, "w") as fout:
        subprocess.run(
            [str(cfg.molscat_exe)],
            stdin=fin, stdout=fout, stderr=subprocess.STDOUT,
            cwd=cfg.work_dir, timeout=timeout,
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
        # GAMMA (ширина Брейта-Вигнера, ПОЛНАЯ -- это НЕ то же самое, что
        # Delta/Gamma_inel из scattering-length методов, единицы совпадают
        # (Гс), но физический смысл другой -- не путать при сведении в
        # общую таблицу) и безразмерная EPSUM_BG (фоновая сумма фаз / pi).
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
                 "Gamma_inel_G", "Gamma_BW_G", "EPSUM_bg_over_pi"):
        result.setdefault(_key, None)

    return result


def characterize_candidate(cfg: Config, candidate: dict, run_id: str) -> dict:
    """Полный цикл характеризации кандидата.

    Стратегия (см. комментарий над регэкспами выше про асимметрию IPRINT
    между IFCONV=2 и IFCONV=3 в реальном коде locpol.f):
      1. IFCONV=1. Если сходится чисто упруго -- готово.
      2. Если LOCPOL просит эскалацию -- пробуем СРАЗУ IFCONV=3 (даёт полный
         набор параметров одним прогоном при IPRINT=6).
      3. NOPEN CHANGED / OSCILLATING / NOT_CONVERGED -- сужаем бракетирующее
         окно вдвое (до 3 раз) и повторяем ТОТ ЖЕ уровень IFCONV.
      4. Если IFCONV=3 так и не сошёлся после всех сужений -- пробуем
         IFCONV=4 (Брейт-Вигнер по сумме собственных фаз) на ИСХОДНОМ окне
         кандидата: это принципиально другая, более устойчивая величина
         (см. переписку про 440 Гс), может сойтись там, где длина рассеяния
         не смогла. ВАЖНО: лимит MXLOC=20 итераций общий для всех IFCONV
         (жёстко зашит в sizes_module.f, не настраивается через &INPUT) --
         смена метода не даёт "больше попыток", только другую, обычно более
         стабильную величину для той же 3-точечной экстраполяции.
      5. Если и это не сошлось -- возвращаем 'failed', но сохраняем ЛУЧШУЮ
         предварительную (не до конца сошедшуюся) оценку B0 из последней
         итерации любого из уровней -- это по-прежнему полезная информация
         (см. 440 Гс: B_RES не сошёлся строго, но стабильно лежал в полосе
         441-444 Гс).

    Все .input/.out этого кандидата складываются в ОТДЕЛЬНУЮ подпапку
    cfg.output_dir/<run_id>/ (а не вперемешку с другими кандидатами) --
    у одного резонанса может быть 4-8 файлов (эскалации + сужения окна),
    подпапка на кандидата держит это читаемым."""
    cand_dir = cfg.output_dir / run_id
    cand_dir.mkdir(parents=True, exist_ok=True)

    fmin, fmax = candidate["field_lo"], candidate["field_hi"]
    ichan = candidate["channel"]  # каждый кандидат несёт свой правильный ICHAN
                                   # (см. detect_candidates -- разные резонансы
                                   # могут сидеть в разных фиксированных каналах)
    ifconv = 1
    narrow_attempts = 0
    escalated = False
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

        if result["status"] == "escalate" and not escalated:
            escalated = True
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

        if ifconv == 3 and not tried_4:
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
                            label: str, tol_cm1: float = 1.0e-8) -> dict:
    """Гоняет ОДНОТОЧЕЧНЫЙ диагностический прогон (FLDMIN=FLDMAX=field) с
    заданным MONQN (через cfg.monqn) и находит номер канала (строку в
    таблице THRESHOLDS), чья энергия совпадает с EREF, посчитанным
    THRSH9 по этим F,mF. ВАЖНО: как обсуждали -- при B=0 разные (F,mF)
    комбинации могут быть точно вырождены (несколько строк с одинаковой
    энергией), поэтому по умолчанию используйте НЕНУЛЕВОЕ поле (см.
    field). Возвращает:
      status='ok': {channel, L, energy_cm1, diff_cm1, ...}
      status='ambiguous': несколько строк совпали в пределах tol_cm1 --
        обычно значит, что field слишком близко к точке вырождения;
        candidates перечисляет все совпавшие строки.
      status='not_found': ни одна строка не совпала -- см. all_thresholds
        для ручной диагностики (например, если MONQN задан некорректно,
        или IPRINT<6 в шаблоне, или ANSA=0 (INUCA=0) -- любой из этих
        случаев даст пустой/бессмысленный список)."""
    if not monqn:
        raise ValueError("find_channel_for_monqn: monqn не задан")
    cfg_local = dataclasses.replace(cfg, monqn=monqn)
    inp = cfg.output_dir / f"{label}.input"
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
                "reason": f"{len(matches)} строк с той же энергией -- "
                "похоже на вырождение (часто бывает у B, близких к 0 "
                "или к точке пересечения уровней); попробуйте другое "
                "поле.", "out_file": str(out)}
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
    coarse_out, coarse_elapsed = run_molscat(cfg, coarse_input, timeout=3600)
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
    p_full.add_argument("--tensor-terms", default="k0i1,k2i1",
                         help="Какие (k,i)-члены V(R) реально включены в этот "
                              "прогон (см. Table 1, Tiesinga et al. NJP 23, "
                              "085007). Через запятую, напр. 'k0i1,k2i1,k2i2' "
                              "чтобы отметить прогон с добавленным магнитным "
                              "диполь-диполем (k=2,i=2). Это ТОЛЬКО метаданные "
                              "для CSV -- какая физика реально считается, "
                              "определяется компиляцией base9-Tm_Tm_AI.f, "
                              "флаг это не переключает.")
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
    p_find.add_argument("--tol-cm1", type=float, default=1.0e-8, dest="tol_cm1",
                         help="Допуск (см⁻¹) при сравнении EREF с "
                              "энергиями из таблицы THRESHOLDS.")

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
            ibfix=args.ibfix, jstep=args.jstep,
        )
        res = find_channel_for_monqn(cfg, args.monqn, args.field, args.label,
                                      tol_cm1=args.tol_cm1)
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
        )
        meta = RunMetadata(
            tensor_terms=args.tensor_terms,
            hyperfine_zeeman=not args.no_hyperfine_zeeman,
            lmax=args.lmax, label=args.label,
        )
        l_filter = None if args.l_filter < 0 else args.l_filter
        run_full_pipeline(cfg, meta, args.fmin, args.fmax, args.dfield,
                           csv_out, l_filter=l_filter,
                           min_separation_G=args.min_separation_g,
                           window_mult=args.window_mult)


if __name__ == "__main__":
    main()
