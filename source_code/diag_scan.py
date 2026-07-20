#!/usr/bin/env python3
"""
diag_scan.py
============
ЛЁГКАЯ диагностика: один грубый скан (без детекции/характеризации) в узком
окне вокруг известного B0, с мелким шагом по полю. Нужен, чтобы быстро
(один вызов MOLSCAT) увидеть, сколько полюсов реально сидит в окрестности
кандидата при данном LMAX -- перед тем как гонять дорогой (минуты) цикл
characterize с автоэскалацией IFCONV.

Переиспользует generate_coarse_scan_input/run_molscat/parse_coarse_scan из
tm_resonance_pipeline.py -- положи этот файл РЯДОМ с tm_resonance_pipeline.py
(в той же source_code/), иначе импорт не найдёт модуль.

Пример:
    python diag_scan.py --exe ./bin/molscat-Tm2.exe \
        --work-dir /c/fortran_projects/molscat/source_code \
        --output-dir /c/fortran_projects/molscat/molscat_runs \
        --fmin 440 --fmax 460 --dfield 0.5 --lmax 6 --ichan 8 --label diag_B447_lmax6
"""

from __future__ import annotations

import argparse
from pathlib import Path

from tm_resonance_pipeline import (
    Config, generate_coarse_scan_input, run_molscat, parse_coarse_scan,
    _fmt_duration,
)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, type=Path)
    ap.add_argument("--work-dir", required=True, type=Path)
    ap.add_argument("--output-dir", type=Path, default=None)
    ap.add_argument("--fmin", type=float, required=True)
    ap.add_argument("--fmax", type=float, required=True)
    ap.add_argument("--dfield", type=float, default=0.5,
                     help="Шаг по полю, Гс (по умолчанию 0.5 -- мелко, чтобы "
                          "не пропустить узкий резонанс вроде Delta=0.18 Гс)")
    ap.add_argument("--lmax", type=int, default=4)
    ap.add_argument("--l", type=int, default=0, dest="l_value",
                     help="Оставить в выводе только строки с этим L "
                          "(по умолчанию 0 -- s-волна). ВНИМАНИЕ: номер CHAN "
                          "не сохраняется между разными LMAX (базис больше -> "
                          "больше каналов -> другая нумерация), поэтому здесь "
                          "НЕТ фильтра по конкретному ICHAN -- печатаются ВСЕ "
                          "каналы с этим L, дальше нужно смотреть глазами, "
                          "какой из них по масштабу re_a похож на нужный фон.")
    ap.add_argument("--expect-abg", type=float, default=None,
                     help="Ожидаемый фон (Re(A), ангстрем) из LMAX=4 прогона "
                          "-- если задан, скрипт подсветит канал(ы), чьи "
                          "значения ближе всего по порядку величины.")
    ap.add_argument("--label", default="diag_scan")
    ap.add_argument("--csv-out", type=Path, default=None,
                     help="Куда сохранить распарсенную таблицу (по умолчанию "
                          "<output-dir>/<label>_diag.csv)")
    args = ap.parse_args()

    args.work_dir.mkdir(parents=True, exist_ok=True)
    output_dir = args.output_dir or (args.work_dir / "molscat_runs")
    output_dir.mkdir(parents=True, exist_ok=True)

    # ВАЖНО: ICHAN для грубого скана НЕ трогаем -- оставляем дефолт
    # cfg.ichan_guess=1, ровно как в run_full_pipeline (см. Config/
    # generate_coarse_scan_input). Это референсный канал для расчёта
    # энергий, а не фильтр -- если сюда подставить номер канала-кандидата
    # (как делала предыдущая версия скрипта), сдвигается порог AWVMAX и
    # меняется сам набор каналов, печатаемых в таблице "low kinetic
    # energy" -- отсюда и пустой результат в прошлый раз.
    cfg = Config(
        molscat_exe=args.exe, work_dir=args.work_dir, output_dir=output_dir,
        coarse_template=Path(""), lmax=args.lmax,
    )

    inp = output_dir / f"{args.label}.input"
    generate_coarse_scan_input(
        cfg, fmin=args.fmin, fmax=args.fmax, dfield=args.dfield,
        label=args.label, out_path=inp,
    )

    print(f"Запускаю MOLSCAT: окно {args.fmin}-{args.fmax} Гс, "
          f"шаг {args.dfield} Гс, LMAX={args.lmax}, ICHAN=1 (референс, "
          f"как в грубом скане Этапа 0)...")
    out_path, elapsed = run_molscat(cfg, inp)
    print(f"Готово за {_fmt_duration(elapsed)} -> {out_path}")

    # Парсим ВСЕ каналы -- номер CHAN не стабилен между разными LMAX,
    # фильтровать по конкретному числу здесь бессмысленно.
    df = parse_coarse_scan(out_path, target_ichan=None)
    if df.empty:
        print("[ПУСТО] Ничего не распарсено -- см. предупреждение выше "
              "(вероятно, IPRINT/энергия не совпадают с ожидаемым).")
        return

    df_l = df[df.L == args.l_value].reset_index(drop=True)
    if df_l.empty:
        print(f"[ПУСТО после фильтра L={args.l_value}] Есть строки для L="
              f"{sorted(df.L.unique())}.")
        df_l = df

    csv_out = args.csv_out or (output_dir / f"{args.label}_diag.csv")
    df_l.to_csv(csv_out, index=False)
    print(f"\nСохранено: {csv_out}\n")

    channels = sorted(df_l.channel.unique())
    print(f"Каналов с L={args.l_value} в этом окне: {channels}\n")

    for ch in channels:
        sub = df_l[df_l.channel == ch].sort_values("field_G")
        re_range = (sub.re_a.min(), sub.re_a.max())
        tag = ""
        if args.expect_abg is not None:
            typical = sub.re_a.abs().median()
            if 0.3 * abs(args.expect_abg) < typical < 3 * abs(args.expect_abg):
                tag = "  <-- по масштабу похож на ожидаемый фон"
        print(f"--- CHAN={ch} (L={args.l_value}), Re(A) диапазон "
              f"[{re_range[0]:.3f}, {re_range[1]:.3f}] Å{tag} ---")
        print(f"{'field_G':>10} {'re_a':>14} {'im_a':>14}   note")
        prev_re = None
        for _, row in sub.iterrows():
            note = ""
            if prev_re is not None:
                jump = abs(row.re_a - prev_re)
                if jump > 3 * max(1.0, abs(prev_re)):
                    note = "<-- резкий скачок Re(A) (возможен полюс рядом)"
            print(f"{row.field_G:10.3f} {row.re_a:14.4f} {row.im_a:14.6f}   {note}")
            prev_re = row.re_a
        print()

    print("Смотри на Re(A) в КАЖДОМ канале отдельно: один плавный переход "
          "через полюс (- -> +-inf -> -) = один резонанс в этом канале; "
          "если у канала, похожего по масштабу на нужный фон, ДВА отдельных "
          "скачка в этом окне -- вот тогда это два близких резонанса.")


if __name__ == "__main__":
    main()
